logit.inv <- function(x) {1/(1+exp(-x))}

estrange <- function(x) {
  rge <- quantile(x,c(.025,.975))
  mn <- x[1]
  yy <- x-x[1]
  ecdf1 <- ecdf(abs(x-mn))
  return(list(rge=rge,mn=mn,ecdf1=ecdf1))
}

predvals <- function(x, ps) {
  lgt <- log(ps/(1-ps))
  lvals <- ifelse(lgt>=x[["rge"]][1] & lgt<=x[["rge"]][2],
                  1-(x[["ecdf1"]](abs(lgt-x[["mn"]]))),
                  0)
  return(lvals)
}

sep <- function(x, ncond, cond) {
  logit.inv <- function(y) {return(1/(1+exp(-y)))}
  sds <- apply(x, 2, sd)
  mnn <- logit.inv(x[1,]-1.96*sds)
  mxx <- logit.inv(x[1,]+1.96*sds)
  return(data.frame(t=as.numeric(colnames(x))+.34,
                    mean=logit.inv(x[1,]),
                    min=mnn,max=mxx,Group=ncond,
                    Comp=cond,
                    Condition=paste(ncond,cond,sep="-")))
}

sep2 <- function(x, ncond, cond) {
  sds <- apply(x, 2, sd)
  mnn <- x[1,]-1.96*sds
  mxx <- x[1,]+1.96*sds
  return(data.frame(t=as.numeric(colnames(x))+.34,
                    mean=x[1,],
                    min=mnn,max=mxx,Group=ncond,
                    Comp=cond,
                    Condition=paste(ncond,cond,sep="-")))
}

calcpred <- function(x, ncond, cond) {
  ff <- apply(x,2,estrange)
  gg <- lapply(ff, predvals, ps=seq(.05,.95,.01))
  hh <- melt(gg)
  hh$x <- as.numeric(hh$L1)+.34
  hh$y <- seq(.05,.95,.01)
  hh$ncond <- ncond
  hh$cond <- cond
  hh$condition <- paste(ncond, cond, sep="-")
  return(subset(hh, value>0))
}
library(ggplot2)
load(file="plottmp.RData")

ff <- sep(NNS.c.lo, ncond="Non-native",cond="control")
gg <- sep(NNS.e.lo, ncond="Non-native",cond="competitor")
hh <- sep(NS.c.lo, ncond="Native",cond="control")
ii <- sep(NS.e.lo, ncond="Native",cond="competitor")

allplot <- rbind(ff,gg,hh,ii)
allplot <- allplot[order(allplot$Group, allplot$Comp, allplot$t),]

NNS.all <- logit.inv(NNS.e.lo)-logit.inv(NNS.c.lo)
NS.all <- logit.inv(NS.e.lo)-logit.inv(NS.c.lo)

jj <- sep2(NNS.all, ncond="Non-native",cond=1)
kk <- sep2(NS.all, ncond="Native",cond=1)
ap2 <- rbind(kk,jj)

ed3.wbin$Condition <- paste(ifelse(ed3.wbin$ncond=="NNS","Non-native","Native"),
														ifelse(ed3.wbin$cond=="control","control","competitor"),sep="-")
ed3.wbin$t <- ed3.wbin$bin/1000
ed3.wbin$time <- ed3.wbin$t
ed3.wbin$probability <- ed3.wbin$p
ed3.wbin <- subset(ed3.wbin, bin<=2500)
ed3.wbin$Group <- ifelse(ed3.wbin$ncond=="NNS","Non-native","Native")
ed3.wbin$Critical <- ifelse(ed3.wbin$cond=="control","Control","Competitor")
ed3.wbin$`time (seconds)` <- ed3.wbin$t
allplot$`time (seconds)` <- allplot$t
allplot$`probability` <- allplot$mean
allplot$`Condition` <- with(allplot, paste(Group, Comp, sep="-"))
allplot$Condition <- factor(allplot$Condition)
allplot$Critical <- ifelse(allplot$Comp == "control", "Control", "Competitor")

ggplot(allplot, aes(x=`time (seconds)`,
										y=`probability`)) +
	geom_ribbon(aes(fill=Group, group = Condition,
									ymin=min, ymax=max), alpha=.2) +
	geom_line(aes(color=Group, shape=Critical)) +
	geom_point(aes(color=Group, shape=Critical), data=ed3.wbin, size = 2) +
	scale_shape_manual(values=c(1,2,1,2)) +
	theme(legend.position="top") +
	labs(shape = "", fill = "", color = "")
