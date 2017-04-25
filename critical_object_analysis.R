library("dplyr")
library("tidyr")
library("purrr")

objs <- load(file = "pogdata.RData")

fit_once <- function(x) {
	x2 <- x %>%
		mutate(O = N - Y,
					 NC = (ncond == "NS") - mean(ncond == "NS"),
					 C = (cond == "exp") - mean(cond == "exp"))

	junk <- capture.output(mod <- nnet::multinom(cbind(Y, O) ~ NC * C, x2))
	coef(mod)["O", ]
}

## permute which group subject belongs to (native vs nonnative)
permute <- function(x, unit = "SubjID") {
	subj_inf <- x %>% select_(unit, "ncond") %>% distinct() %>%
		mutate(ncond = sample(ncond))
	inner_join(x %>% select(-ncond), subj_inf, unit)
}

## permute experimental versus control condition (by subj)
permute_exp <- function(x) {
	reverse_it <- function(x) {
		x[["cond"]] <- rev(x[["cond"]])
		x
	}
	## use synchronized permutation logic (pesarin/salmaso)
	x_nest <- x %>% nest(-SubjID, -ncond, .key = "dat") 
	x_nest_g <- x_nest %>%
		group_by(ncond)
	stopifnot(length(unique(group_size(x_nest_g))) == 1L)
	cflip <- sample(c(TRUE, FALSE), group_size(x_nest_g)[1], replace = TRUE)
	ff <- x_nest_g %>%
		mutate(flip = sample(cflip),
					 dat = map_if(dat, flip, reverse_it)) %>%
		ungroup() %>%
		unnest()
}

permute_item <- function(x) {
	flip_coin <- function(xx) {
		scheme <- matrix(c(1, 2, 3, 4,
											 3, 4, 1, 2), ncol = 2)
		xx[scheme[, sample(1:2, 1)], "oldcond"] %>%
			rename(ncond = oldcond)
	}
	item_inf <- x %>% select(ItemID, oldcond = ncond, cond) %>%
		nest(-ItemID, .key = "oldcond") %>%
		mutate(newcond = map(oldcond, flip_coin)) %>%
		unnest(old = oldcond, newcond)
	x %>% rename(oldcond = ncond) %>%
		inner_join(item_inf, c("ItemID", "oldcond", "cond")) %>%
		select(-oldcond)
}

## permute experimental versus control condition (by item)
permute_exp_item <- function(x) {
	shuf_tbl <- x %>%
		distinct_("ItemID", "cond") %>%
		rename(cond2 = cond) %>%
		group_by_("ItemID") %>%
		mutate(cond = sample(cond2)) %>%
		ungroup()
	inner_join(shuf_tbl, x, c("ItemID", "cond2" = "cond")) %>%
		select(-cond2)
}


cutoff <- quantile(ri$RT, .975)

## pull out the analysis window
## goes from 200 ms after speech onset until end of trial
## or 'cutoff' point (whichever is earlier)
pog2 <- subset(pogm, (ms >= 200) & (ms <= ifelse(RT > cutoff, cutoff, RT)) )

## convert "ID" (the response category) to 1 or 0
## depending on whether there is a gaze to the critical object
pog2$crit <- ifelse(pog2$ID == "C", 1, 0)

## ## ## analysis by subject

gaze_counts <- pog2 %>%
	group_by(SubjID, ncond, cond) %>%
	summarize(Y = sum(crit), N = n()) %>%
	ungroup()

## p-values for main effect of NS/NNS & native:cond interaction
orig <- fit_once(gaze_counts)
mx <- cbind(orig, replicate(9999L, permute(gaze_counts) %>%
																	 fit_once()))

p_vals <- apply(mx[c(2, 4), ], 1,
								function(x) sum(abs(x) >= abs(x[1]))) / ncol(mx)

## get p-value for main effect of exp/control condition
mx_cond <- cbind(orig, replicate(9999L,
																 permute_exp(gaze_counts) %>%
																 fit_once()))

p_val_cond <- sum(abs(mx_cond[3, ]) >= abs(mx_cond[3, 1])) /
	length(mx_cond[3, ])

p1 <- c(p_vals[1], C = p_val_cond, p_vals[2])

## ## ## analysis by item

gaze_counts_item <- pog2 %>%
	group_by(ItemID, ncond, cond) %>%
	summarize(Y = sum(crit), N = n()) %>%
	ungroup()

## p-values for main effect of NS/NNS & native:cond interaction
orig_item <- fit_once(gaze_counts_item)

mx_item <- cbind(orig_item,
								 replicate(9999L, permute_item(gaze_counts_item) %>%
																	fit_once()))
p_vals_item <- apply(mx_item[c(2, 4), ], 1,
										 function(x) sum(abs(x) >= abs(x[1]))) /
	ncol(mx_item)

## get p-value for main effect of exp/control condition
mx_cond_item <-
	cbind(orig, replicate(9999L,
												permute_exp_item(gaze_counts_item) %>%
												fit_once()))

p_val_cond_i <- sum(abs(mx_cond_item[3, ]) >=
										abs(mx_cond_item[3, 1])) /
	length(mx_cond_item[3, ])

p2 <- c(p_vals_item[1], C = p_val_cond_i, p_vals_item[2])

results <- rbind(param_est = orig,
								 param_nns = c(orig[1] + -.5 * orig[2], NA,
															 orig[3] + -.5 * orig[4], NA),
								 param_ns = c(orig[1] + .5 * orig[2], NA,
															orig[3] + .5 * orig[4], NA),
								 p1 = c(NA, p1),
								 p2 = c(NA, p2))

## NC: effect of native vs nonnative speaker
## C:  effect of competitor vs noncompetitor
## 
## exp(-results["param_ns", "C"]) = 5.81; odds ratio Native
## exp(-results["param_nns", "C"]) = 4.19; odds ratio non-native

saveRDS(results, "results_competitor_analysis.rds")
