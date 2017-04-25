library(abind)
library(plyr)

# x is a vector of frame data
crBaselineCorrect <- function(x, frameEnd) {
		blcorrect <- function(x, f1) {
				ff <- x-mean(x[1:f1])
				ff[-(1:f1)]
		}
		alldims <- 1:length(dim(x))
		dims <- alldims[-2]
		newmatrix <- apply(x, dims, blcorrect, f1=frameEnd)
		aperm(newmatrix, c(2,1,setdiff(alldims,1:2)))
}

# calculate signed negative log of p for each parameter and each frame
# ax : array of permutation values, first index is run, second index is frame
# dimensions 3...length(dim(pax)) hold the parameters
crTestStats <- function(ax, baselineEnd=NULL) {
		# snlogp() gives the signed negative log of p, where
		# ix : which run to calculate clusters on (1=original)
		snlogp <- function(x, ix=1) {
				# calculate proportion greater than or equal to original (ix)
				-log(sum(abs(x)>=abs(x[ix]), na.rm=TRUE)/length(x))*sign(x[ix])
		}
		# TODO: baseline correction
		if (!is.null(baselineEnd)) {
				ax <- crBaselineCorrect(ax, baselineEnd)
		} else {}
		last_dim <- length(dim(ax))  # frame dimension
		nruns <- dim(ax)[1] # number of runs
		# treat each run as 'original', calculate p-values
		ff <- laply(1:nruns, function(x) {
				# for each parameter/frame combination
				apply(ax, 2:last_dim, snlogp, ix=x)
		})
}

# given a vector x, pull out all significant clusters
# and return either the max cluster mass statistic, or
# a dataframe listing all the clusters
crGetClusters <- function(x, cutoff=.05, maxcmsOnly=TRUE) {
		# filter out all the nonsignificant clusters
		vec <- (abs(x)>=(-log(cutoff)))*x
		# next line is just for testing
		# vec <- c(rep(5,4),rep(0,7),4,rep(0,3),rep(-5,2),rep(4.1,1))
		runs <- rle(sign(vec))
		clust_index <- (1:length(runs$lengths))[runs$values!=0]
		# calculate onset (t0), offset (t1) and cluster mass stat (cms)
		# for all clusters
		result <- ldply(clust_index, function(lx) {
				# get start of cluster (element index)
				if (lx>1) {
						t0 <- sum(runs$lengths[1:(lx-1)])+1
				} else {
						t0 <- 1
				}
				# end of cluster (element index)
				t1 <- t0+runs$lengths[lx]-1
				data.frame(t0=t0, t1=t1, cms=abs(sum(vec[t0:t1])))
		})
		if (maxcmsOnly) {
				if (nrow(result)>0) {
						result <- max(result$cms)
				} else {
						result <- 0
				}
		} else {}
		return(result)
}

# get cluster mass statistics for each run
# ax : array of cluster test statistics
#      dimension 1 is run, 2nd dimension is frame number
#      remaining dims are the different parameters
crAnalyze <- function(ax, cutoff=.05) {
		last_dim <- length(dim(ax))
		nruns <- dim(ax)[1]
		adims <- (1:length(dim(ax)))[-2] # all dimensions, except frame
		apply(ax, adims, crGetClusters, cutoff=cutoff)
}

# compare list of clusters (for each parameter)
# to distribution (for each parameter)
crResults <- function(clustlist, maxcms) {
		lapply(1:length(clustlist), function(lx) {
				if (nrow(clustlist[[lx]])==0) { # no clusters detected
						result <- NULL
				} else {
						nge <- sapply(clustlist[[lx]]$cms, function(x) {
								sum(maxcms[,lx]>=x)
						})
						result <- cbind(clustlist[[lx]], nge=nge, pval=nge/length(maxcms[,lx]))
				}
				return(result)
		})
}  

# randomize conditions
scramble <- function(dat, newlabels) {
		newdat <- merge(dat, newlabels, by="SubjID")
}

fitOnce <- function(dat) {
		library(plyr)
		library(nnet)
		dat$C <- ifelse(dat$cond=="exp",.5,-.5)
		dat$S <- ifelse(dat$ncond2=="NS",.5,-.5)
		res <- daply(dat, .(bin), function(x) {
				capture.output(ff <- coef(multinom(cbind(crit,ncrit)~C*S, x, maxit=1000)))
				ff
		})
		res
}

doOnce <- function(dat, subj) {
		subj$ncond2 <- sample(subj$ncond2) 
		fitOnce(scramble(dat, subj))
}

doOnceComp <- function(dat, subj) {
		subj2 <- ddply(subj, .(SubjID), function(x) {
				x$cond2 <- sample(x$cond)
				x
		})
		ff <- merge(dat, subj2, by=c("SubjID","cond"))
		ff <- ff[,setdiff(colnames(ff),"cond")]
		colnames(ff) <- sub("cond2","cond",colnames(ff))
		colnames(ff) <- sub("ncond","ncond2",colnames(ff))
		fitOnce(ff)
}

doOnceItem <- function(dat, item) {
		item2 <- ddply(item, .(ItemID), function(x) {
				x$ncond2 <- sample(x$ncond)
				x
		})
		dat.item <- merge(dat, item2, by=c("ItemID","ncond"))
		fitOnce(dat.item)
}

load(file="pogdata.RData", verbose=TRUE)
# pogm contains eye data
# ri contains information about each trial

# truncate trials at top 2.5% of the RT distribution
# cutoff <- quantile(ri$RT, .975)
cutoff <- 2500

## number of monte carlo simulations
nmc <- 9999L

# pull out the analysis window
# goes from 200 ms after speech onset until end of trial
# or 'cutoff' point (whichever is earlier)
pog2 <- subset(pogm, (ms >= 200) & (ms <= ifelse(RT>cutoff,cutoff,RT)) )
pog2$SubjID <- as.numeric(as.character(pog2$SubjID))
# convert "ID" (the response category) to 1 or 0
# depending on whether there is a gaze to the critical object
pog2$crit <- ifelse(pog2$ID=="C",1,0)

pog2$bin <- floor( ((pog2$Frame-101)+12)/24 )*48+200
pog3.agg <- aggregate(crit ~ bin + cond + ncond + SubjID, pog2, sum)
pog3.n <- aggregate(crit ~ bin + cond + ncond + SubjID, pog2, length)
colnames(pog3.n) <- sub("crit","n",colnames(pog3.n))
pog3 <- merge(pog3.agg, pog3.n)
pog3$ncrit <- pog3$n-pog3$crit
pog3 <- pog3[order(pog3$SubjID, pog3$cond, pog3$bin),]

subj <- subset(as.data.frame(xtabs(~ncond+SubjID, pog3)), Freq>0)
subj <- subj[order(subj$SubjID, subj$ncond),c("SubjID","ncond")]
colnames(subj) <- sub("ncond", "ncond2", colnames(subj))

bins <- unique(pog3$bin)

# analysis of Speaker
pog3 <- subset(pog3, bin<=1016)
runs <- raply(nmc, doOnce(pog3, subj), .progress="text")

orig.bySubj <- array(fitOnce(scramble(pog3, subj)), c(1,18,4))
allruns <- abind(orig.bySubj, runs, along=1)

apply(allruns, 2:3, function(x) {sum(abs(x)>=abs(x[1]))/length(x)})

crtest.bySubj <- crTestStats(allruns[,,-c(1,2)])
crstat.bySubj <- crAnalyze(crtest.bySubj)
crOrig.bySubj <- apply(crtest.bySubj[1,,], 2, crGetClusters, maxcmsOnly=FALSE)
res <- crResults(crOrig.bySubj, crstat.bySubj)
results.bySubj <- rbind("Speaker" = res[[1]],
												"Speaker:Condition" = res[[2]])

# create data for item analysis
pog3.item.agg <- aggregate(crit ~ bin + cond + ncond + ItemID, pog2, sum)
pog3.item.n <- aggregate(crit ~ bin + cond + ncond + ItemID, pog2, length)
colnames(pog3.item.n) <- sub("crit","n",colnames(pog3.item.n))
pog3.item <- merge(pog3.item.agg, pog3.item.n)
pog3.item$ncrit <- pog3.item$n-pog3.item$crit
pog3.item <- pog3.item[order(pog3.item$ItemID, pog3.item$cond, pog3.item$bin),]

item <- subset(as.data.frame(xtabs(~ncond+ItemID, pog3.item)), Freq>0)
item <- item[order(item$ItemID, item$ncond),c("ItemID","ncond")]
item2 <- item
item2$ncond2 <- item2$ncond

# analysis of Speaker by Item
pog3.byItem <- subset(pog3.item, bin<=1016)
runs.byItem <- raply(nmc, doOnceItem(pog3.byItem, item), .progress="text")

orig.byItem <- array(fitOnce(merge(pog3.byItem, item2, by=c("ItemID","ncond"))), c(1,18,4))
allruns.byItem <- abind(orig.byItem, runs.byItem, along=1)

apply(allruns.byItem, 2:3, function(x) {sum(abs(x)>=abs(x[1]))/length(x)})

crtest.byItem <- crTestStats(allruns.byItem[,,-c(1,2)])
crstat.byItem <- crAnalyze(crtest.byItem)
crOrig.byItem <- apply(crtest.byItem[1,,], 2, crGetClusters, maxcmsOnly=FALSE)
res <- crResults(crOrig.byItem, crstat.byItem)
results.byItem <- rbind("Speaker" = res[[1]],
												"Speaker:Condition" = res[[2]])

# OK, now repeat, doing the analysis of competitor
subj.comp <- subset(as.data.frame(xtabs(~cond+SubjID, pog3)), Freq>0)
subj.comp <- subj.comp[order(subj.comp$SubjID, subj.comp$cond),c("SubjID","cond")]

runs.comp.bySubj <- raply(nmc, doOnceComp(pog3, subj.comp), .progress="text")
pog.tmp <- pog3
colnames(pog.tmp) <- sub("^ncond$","ncond2",colnames(pog.tmp))
colnames(pog.tmp) <- sub("^cond$","cond2",colnames(pog.tmp))

orig.comp.bySubj <- array(fitOnce(pog.tmp), c(1,18,4))
allruns.comp.bySubj <- abind(orig.comp.bySubj, runs.comp.bySubj, along=1)

apply(allruns.comp.bySubj, 2:3, function(x) {sum(abs(x)>=abs(x[1]))/length(x)})

crtest.comp.bySubj <- crTestStats(allruns.comp.bySubj[,,-c(1,3)])
crstat.comp.bySubj <- crAnalyze(crtest.comp.bySubj)
crOrig.comp.bySubj <- apply(crtest.comp.bySubj[1,,], 2, crGetClusters, maxcmsOnly=FALSE)
res <- crResults(crOrig.comp.bySubj, crstat.comp.bySubj)
results.comp.bySubj <- rbind("Condition" = res[[1]])

# Competition, by Item
pog.itm.tmp <- pog3.byItem
colnames(pog.itm.tmp) <- sub("ItemID","SubjID",colnames(pog.itm.tmp))
item.comp <- subset(as.data.frame(xtabs(~cond+SubjID, pog.itm.tmp)), Freq>0)
item.comp <- item.comp[order(item.comp$SubjID, item.comp$cond),c("SubjID","cond")]

runs.comp.byItem <- raply(nmc, doOnceComp(pog.itm.tmp, item.comp), .progress="text")
pog.tmp <- pog.itm.tmp
colnames(pog.tmp) <- sub("^ncond$","ncond2",colnames(pog.tmp))
colnames(pog.tmp) <- sub("^cond$","cond2",colnames(pog.tmp))

orig.comp.byItem <- array(fitOnce(pog.tmp), c(1,18,4))
allruns.comp.byItem <- abind(orig.comp.byItem, runs.comp.byItem, along=1)

apply(allruns.comp.byItem, 2:3, function(x) {sum(abs(x)>=abs(x[1]))/length(x)})

crtest.comp.byItem <- crTestStats(allruns.comp.byItem[,,-c(1,3)])
crstat.comp.byItem <- crAnalyze(crtest.comp.byItem)
crOrig.comp.byItem <- apply(crtest.comp.byItem[1,,], 2, crGetClusters, maxcmsOnly=FALSE)
res <- crResults(crOrig.comp.byItem, crstat.comp.byItem)
results.comp.byItem <- rbind("Condition" = res[[1]])

all_results <- cbind(Effect = c(rownames(results.bySubj),
								 rownames(results.comp.bySubj),
								 rownames(results.byItem),
								 rownames(results.comp.byItem)),
			Unit = rep(c("Subject", "Item"),
								 c(sum(nrow(results.bySubj), nrow(results.comp.bySubj)),
									 sum(nrow(results.byItem), nrow(results.comp.byItem)))),
			rbind(results.bySubj, results.comp.bySubj, results.byItem,
						results.comp.byItem))

rownames(all_results) <- NULL
all_results$b0 <- bins[all_results$t0]
all_results$b1 <- bins[all_results$t1]
all_results <- all_results[order(all_results$Effect, all_results$Unit), ]

print(all_results)

saveRDS(all_results, "results_cluster_randomization.rds")
