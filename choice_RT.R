library(lme4)

rih <- readRDS("choice_data.rds")

xtabs(~choice+cond+ncond, rih)

rih$err <- ifelse(rih$choice!="TARGET",1,0)
rih.comp <- subset(rih, cond=="exp")
rih.comp$spkr <- scale(ifelse(rih.comp$ncond=="NNS",0,1), scale=F)

m1 <- glmer(err ~ spkr + (1 | SubjID) + (1 + spkr | ItemID),
						data=rih.comp,
						family=binomial)

m2 <- glmer(err ~ (1 | SubjID) + (1 + spkr | ItemID),
						data=rih.comp,
						family=binomial)

anova(m1, m2)

rtcut <- quantile(rih$RT, .975)
rih$RTt <- ifelse(rih$RT>=rtcut, rtcut, rih$RT)

rih$C <- scale(ifelse(rih$cond=="control",0,1), scale=FALSE)
rih$S <- scale(ifelse(rih$ncond=="NNS",1,0),scale=FALSE)

rt.m1 <- lmer(RTt ~ C*S + (1 + C | SubjID) + (1 + C*S | ItemID),
							data=rih, REML=FALSE)
rt.noCS <- lmer(RTt ~ C*S-C:S + (1 + C | SubjID) + (1 + C*S | ItemID),
								data=rih, REML=FALSE)
rt.noC <- lmer(RTt ~ C*S-C + (1 + C | SubjID) + (1 + C*S | ItemID),
								data=rih, REML=FALSE)
rt.noS <- lmer(RTt ~ C*S-S + (1 + C | SubjID) + (1 + C*S | ItemID),
								data=rih, REML=FALSE)

anova(rt.m1, rt.noCS)
anova(rt.m1, rt.noC)
anova(rt.m1, rt.noS)

with(rih, aggregate(list(RT=RTt), list(cond=cond, ncond=ncond), mean))
with(rih, aggregate(list(RT=RTt), list(cond=cond), mean))
with(rih, aggregate(list(RT=RTt), list(ncond=ncond), mean))

with(rih, aggregate(list(RT=RTt), list(cond=cond, ncond=ncond), sd))
with(rih, aggregate(list(RT=RTt), list(cond=cond), sd))
with(rih, aggregate(list(RT=RTt), list(ncond=ncond), sd))
