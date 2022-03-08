library(vegan)
#install.packages("remotes")
remotes::install_github("jfq3/QsRutils")
library(QsRutils)
library(tidyr)
library(dplyr)
library(ecodist)
library(ape)
library(lme4)
library(ggplot2)
library(multcomp)
library(emmeans)
library(nlme)
library(plyr)
library(SingleCaseES)
library(WebPower)



# Study metadata
meta <- read.csv("Raw-data/study-sequences/Gao/SraRunTable_updated.csv")
# Study microbial diversity data
divdat <- read.csv("Processed-data/Gao-diversity.csv")


# combine data into one dataset
dat <- merge(meta, divdat, by.x="Experiment", by.y="ID")

dat$Cropping.system <- as.factor(dat$Cropping.system)
dat$Cropping.system <- factor(dat$Cropping.system, levels(dat$Cropping.system)[c(4,2,3,1)])
dat$Replicate <- as.factor(dat$Replicate.y)


dat[,c("Sample.Name", "Cropping.system")]







colorvec <- c("burlywood4", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4")



###### BACTERIA ###### 

################  MICROBIAL DIVERSITY RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system")) 
options("scipen"=100, "digits"=4)



## inverse Simpson's diversity ## 
dat$invsimpson.div <- as.numeric(dat$invsimpson.div)
hist((dat$invsimpson.div))
# mixed effects model
fit <- lme(invsimpson.div ~ Cropping.system, random = ~1|Replicate, data = dat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat[which(is.na(dat$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system)%>% 
  summarise_at(vars(invsimpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$invsimpson.div <- ssp


# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(dat, aes(x=Cropping.system, y=invsimpson.div, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Fungal Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat$invsimpson.div)) + c(0,30))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1,2,3,4), y = rep(225, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Gao_agroforestry.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Gao_diversity_anova.csv")



# From this figure it looks like there could be an effect of AF that increases with duration.
# Based on treatment averages, the 5, 9, and 14-year AF treatments cause the following percent increases in diversity:
af5 <- ((testlet$emmean[2]-testlet$emmean[1])/testlet$emmean[1])*100 
af9 <- ((testlet$emmean[3]-testlet$emmean[1])/testlet$emmean[1])*100 
af14 <- ((testlet$emmean[4]-testlet$emmean[1])/testlet$emmean[1])*100 


# calculate group st dev
resid <- residuals(fit)
resid = cbind(divdat[which(is.na(divdat$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system)%>% 
  summarise_at(vars(invsimpson.div),funs(stdev=sd(., na.rm=T)))

# calculate effect sizes
af5_es <-((testlet$emmean[2]-testlet$emmean[1])/mean(temp$stdev[which(temp$Cropping.system %in% c("Wheat", "AF5"))]))
af9_es <-((testlet$emmean[3]-testlet$emmean[1])/mean(temp$stdev[which(temp$Cropping.system %in% c("Wheat", "AF9"))]))
af14_es <-((testlet$emmean[4]-testlet$emmean[1])/mean(temp$stdev[which(temp$Cropping.system %in% c("Wheat", "AF14"))]))



### power-sample size curves
ptest5 <- wp.t(n1=seq(2,200,1), d = af5_es, alpha = 0.05, alternative="two.sided")
ptest9 <- wp.t(n1=seq(2,200,1), d = af9_es, alpha = 0.05, alternative="two.sided")
ptest14 <- wp.t(n1=seq(2,200,1), d = af14_es, alpha = 0.05, alternative="two.sided")

pdat <- data.frame(treatment = levels(dat$Cropping.system)[c(2:4)],
                   percent.change = round(c(af5, af9, af14),2), 
                   es = round(c(af5_es, af9_es, af14_es),2),
                   n80 = c(ptest5$n[min(which(ptest5$power>0.8))], ptest9$n[min(which(ptest9$power>0.8))], ptest14$n[min(which(ptest14$power>0.8))]))
write.csv(pdat, paste0("Model-output/power-tests/Gao-invshannon.div.csv"))


pdf("Figures/effect-plots/*3_power-test_Gao.pdf", height=4, width=5)
par(mar=c(4,4,0.5,0.5))
plot(ptest5,type='b', las=1, bg=alpha("darkseagreen2",0.8), ylim=c(-0.05,1.01), xlim=c(-1,ptest5$n[min(which(ptest5$power>0.8))])+5, pch=21, col=NA)
points(ptest9$power~ptest9$n, type='b', las=1, bg=alpha("darkseagreen3",0.8), pch=21, col=NA)
points(ptest14$power~ptest14$n, type='b', las=1, bg=alpha("darkseagreen4",0.8), pch=21, col=NA)
legend(60, 0.4, legend=c("AF14", "AF9", "AF5"), pt.bg=c("darkseagreen4", "darkseagreen3", "darkseagreen2"), pch=21)
points(ptest5$n[min(which(ptest5$power>0.8))], 0, pch=18, col=alpha("darkseagreen2",0.8), cex=4)
points(ptest9$n[min(which(ptest9$power>0.8))], 0, pch=18, col=alpha("darkseagreen3",0.8), cex=4)
points(ptest14$n[min(which(ptest14$power>0.8))], 0, pch=18, col=alpha("darkseagreen4",0.8), cex=4)
text(pdat$n80, y=0, labels=paste0("n=",pdat$n80), cex=1)
abline(h=0.8, lty=2)
segments(x0=c(5,16,103), x1=c(5,16,103), y0=0.8, y1=0.075, lty=2)
#mtext("C", side=3, line=0.5, cex=1.1, adj = 0)
dev.off()









## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

numcomp <- 3
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "AF5" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$invsimpson.div[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$invsimpson.div[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "AF9" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$invsimpson.div[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$invsimpson.div[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 3
comp <- 3
Lt <- "AF14" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$invsimpson.div[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$invsimpson.div[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Gao_diversity_RR.csv")








################  CARBON RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system")) 
options("scipen"=100, "digits"=4)




## SOC ## 

testlet <- aggregate(SOC~Cropping.system, dat=dat, FUN="mean")
testlet_SOC <- testlet[order(testlet$Cropping.system),] 

p <- ggplot(dat, aes(x=Cropping.system, y=SOC, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic carbon (g kg"^-1, ")"))) +
  lims(y=c(range(na.omit(dat$SOC)) + c(0,0.1))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) #+ geom_text(x = c(1,2,3,4), y = rep(3.3, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Gao_SOC_cropsys.pdf")












## DOC ## 

testlet <- aggregate(DOC~Cropping.system, dat=dat, FUN="mean")
testlet_DOC <- testlet[order(testlet$Cropping.system),] 

p <- ggplot(dat, aes(x=Cropping.system, y=DOC, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Dissolved organic carbon (g kg"^-1, ")"))) +
  lims(y=c(range(na.omit(dat$DOC)) + c(0,0.1))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) #+ geom_text(x = c(1,2,3,4), y = rep(3.3, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Gao_DOC_cropsys.pdf")






## effect sizes for carbon response to crop diversification ##
# use log-response ratio 

# SOC # 
numcomp <- 3
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "AF5" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$SOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$SOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)


# comparison 2
comp <- 2
Lt <- "AF9" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$SOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$SOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 4)


# comparison 3
comp <- 3
Lt <- "AF14" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$SOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$SOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 4)




write.csv(d, "Model-output/effects_diversity/Gao_carbon-SOC_RR.csv")





# DOC # 
numcomp <- 3
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "AF5" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$DOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$DOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)


# comparison 2
comp <- 2
Lt <- "AF9" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$DOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$DOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 4)


# comparison 3
comp <- 3
Lt <- "AF14" # label for diversity treatment
Lc <- "Wheat" # label for diversity control
Yt <- dat$DOC[which(dat$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat$DOC[which(dat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 4)




write.csv(d, "Model-output/effects_diversity/Gao_carbon-DOC_RR.csv")



################  MICROBIAL DIVERSITY ASSOCIATION WITH SOIL CARBON  ################  

# scale x and y variables to get standardized beta coefficients


lmtest_DOC <- lme(scale(DOC) ~ scale(invsimpson.div), random = ~1|Replicate, data = dat, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_DOC)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_DOC_sum <- as.data.frame(summary(lmtest_DOC)$tTable)
dat$predlmDOC[which(is.na(dat$DOC)==F)] <- predict(lmtest_DOC, level=0)


lmtest_SOC <- lme(scale(SOC) ~ scale(invsimpson.div), random = ~1|Replicate, data = dat, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_SOC)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_SOC_sum <- as.data.frame(summary(lmtest_SOC)$tTable)
dat$predlmSOC[which(is.na(dat$SOC)==F)] <- predict(lmtest_SOC, level=0)



df <- data.frame(statistic = c("intercept", "slope", "t", "df", "p"),
                 value_DOC = c(lmtest_DOC_sum$Value[1], lmtest_DOC_sum$Value[2], lmtest_DOC_sum$`t-value`[2], lmtest_DOC_sum$DF[2], lmtest_DOC_sum$`p-value`[2]),
                 value_SOC = c(lmtest_SOC_sum$Value[1], lmtest_SOC_sum$Value[2], lmtest_SOC_sum$`t-value`[2], lmtest_SOC_sum$DF[2], lmtest_SOC_sum$`p-value`[2]))


write.csv(df, "Model-output/effects_diversity/Gao_lmm-carbon.csv")






p <- ggplot(dat, aes(x=invsimpson.div, y=SOC)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="Bacterial Inverse Simpson diversity", y=expression(paste("Soil organic carbon (g kg"^-1, ")")), fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cropping.system)) +
  theme_bw()  +
  scale_fill_manual(values=colorvec) +
  theme(axis.title=element_text(size=10))
p


q <- ggplot(dat, aes(x=invsimpson.div, y=DOC)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="Bacterial Inverse Simpson diversity", y=expression(paste("Dissolved organic carbon (g kg"^-1, ")")), fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cropping.system)) +
  theme_bw()  +
  scale_fill_manual(values=colorvec) +
  theme(axis.title=element_text(size=10))
q



r <- ggpubr::ggarrange(p+ ggpubr::rremove("xlab"),q+ ggpubr::rremove("xlab"),
                       nrow=1, common.legend=TRUE, legend="right")
s <- ggpubr::annotate_figure(r, bottom = ggpubr::text_grob("Bacterial Inverse Simpson diversity"))
ggpubr::ggexport(s, filename="Figures/regressions/Gao_C-invsimpson.div_bac.pdf", height=3.25,width=8)






