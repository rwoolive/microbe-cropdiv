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



# metadata including diversity calculations
meta <- read.csv("Raw-data/study-sequences/Cloutier-2020/CoverCropFungal-CarbonData.csv")
meta$X
divdat <- read.csv("Processed-data/Cloutier-diversity_amf.csv")
divdat$Season_num <- divdat$Season
divdat$Season_num <- as.factor(recode(divdat$Season_num, Spring = "1", Summer = "2"))
divdat$Cover2 <- divdat$Cover
divdat$Cover2 <- as.factor(recode(divdat$Cover2, "3 Spp mix" = "3SppN", "6 Spp mix" = "6Spp"))
divdat$X <- paste0(divdat$Cover2, "-", divdat$Replicate, "-", divdat$Season_num)
divdat$X

# combine data into one dataset
divdat <- merge(meta, divdat, by="X")
divdat$Cover <- as.factor(divdat$Cover)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(5,3:4,6:9,1,2)])
divdat$Season <- as.factor(divdat$Season)
divdat$Replicate <- as.factor(divdat$Replicate)



colorvec <- c("burlywood4", rep("darkolivegreen2",6), "darkolivegreen3", "darkolivegreen4")



################  MICROBIAL DIVERSITY RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cover, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cover",
                                   "Season",
                                   "Cover:Season")) 
options("scipen"=100, "digits"=4)




## inverse Simpson's diversity ## 
divdat$invsimpson.div <- as.numeric(divdat$invsimpson.div)
hist(sqrt(divdat$invsimpson.div))
# mixed effects model
fit <- lme(sqrt(invsimpson.div) ~ Cover*Season, random = ~1|Replicate, data = divdat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(divdat[which(is.na(divdat$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cover, Season)%>% 
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
test1 <- emmeans(fit, ~Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cover),]
p <- ggplot(divdat, aes(x=Cover, y=invsimpson.div, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cover", y=expression(paste("AMF Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(divdat$invsimpson.div)) + c(0,3))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:9), y = rep(13, 9), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=7, filename = "Figures/effect-plots/Cloutier_fun_AMF_cover.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Cloutier_diversity_AMF_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

numcomp <- 8
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Canola" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 2
comp <- 2
Lt <- "Clover" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 3
comp <- 3
Lt <- "Oat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 4
comp <- 4
Lt <- "Pea" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 5
comp <- 5
Lt <- "Radish" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 6
comp <- 6
Lt <- "Rye" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 7
comp <- 7
Lt <- "3 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 8
comp <- 8
Lt <- "6 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)






write.csv(d, "Model-output/effects_diversity/Cloutier_diversity_AMF_RR.csv")








################  CARBON RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cover, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cover",
                                   "Season",
                                   "Cover:Season")) 
options("scipen"=100, "digits"=4)




## POXC ## 
divdat$POXC <- as.numeric(divdat$POXC)
hist((divdat$POXC))
# mixed effects model
fit <- lme((POXC) ~ Cover*Season, random = ~1|Replicate, data = divdat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(divdat[which(is.na(divdat$POXC)==FALSE),],resid)
temp = resid%>%
  group_by(Cover,Season)%>% 
  summarise_at(vars(POXC),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$POXC <- ssp

# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet_POXC <- testlet[order(testlet$Cover),]
p <- ggplot(divdat, aes(x=Cover, y=POXC, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Permanganate oxidizable carbon (mg kg"^-1,")"))) +
  lims(y=c(range(na.omit(divdat$POXC)) + c(0,30))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:9), y = rep(420, 9), aes(label = .group), data = testlet_POXC)
p
ggpubr::ggexport(p, height=4, width=7, filename = "Figures/effect-plots/Cloutier_POXC_cover.pdf")



## OM ## 
divdat$OM <- as.numeric(divdat$OM)
hist((divdat$OM)^3)
# mixed effects model
fit <- lme((OM)^3 ~ Cover*Season, random = ~1|Replicate, data = divdat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(divdat[which(is.na(divdat$OM)==FALSE),],resid)
temp = resid%>%
  group_by(Cover,Season)%>% 
  summarise_at(vars(OM),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$OM <- ssp

# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet_OM <- testlet[order(testlet$Cover),]
p <- ggplot(divdat, aes(x=Cover, y=OM, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic matter (%)"))) +
  lims(y=c(range(na.omit(divdat$OM)) + c(0,0.1))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:9), y = rep(2.75, 9), aes(label = .group), data = testlet_OM)
p
ggpubr::ggexport(p, height=4, width=7, filename = "Figures/effect-plots/Cloutier_OM_cover.pdf")




write.csv(t(mod.effects), "Model-output/effects_diversity/Cloutier_carbon_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 


# OM #
numcomp <- 8
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Canola" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 2
comp <- 2
Lt <- "Clover" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 3
comp <- 3
Lt <- "Oat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 4
comp <- 4
Lt <- "Pea" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 5
comp <- 5
Lt <- "Radish" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 6
comp <- 6
Lt <- "Rye" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 7
comp <- 7
Lt <- "3 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 8
comp <- 8
Lt <- "6 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$OM[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$OM[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)





write.csv(d, "Model-output/effects_diversity/Cloutier_carbon-OM_RR.csv")







# POXC #
numcomp <- 8
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Canola" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 2
comp <- 2
Lt <- "Clover" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 3
comp <- 3
Lt <- "Oat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 4
comp <- 4
Lt <- "Pea" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 5
comp <- 5
Lt <- "Radish" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 6
comp <- 6
Lt <- "Rye" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 7
comp <- 7
Lt <- "3 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 8
comp <- 8
Lt <- "6 Spp mix" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- divdat$POXC[which(divdat$Cover==Lt)] # diversity treatment mean
Yc <- divdat$POXC[which(divdat$Cover==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)





write.csv(d, "Model-output/effects_diversity/Cloutier_carbon-POXC_RR.csv")






################  MICROBIAL DIVERSITY ASSOCIATION WITH SOIL CARBON  ################  
# scale x and y variables to get standardized beta coefficients


lmtest_POXC <- lme(scale(POXC) ~ scale(invsimpson.div), random = ~1|Replicate, data = divdat, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_POXC)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_POXC_sum <- as.data.frame(summary(lmtest_POXC)$tTable)
divdat$predlmPOXC[which(is.na(divdat$POXC)==F)] <- predict(lmtest_POXC, level=0)


lmtest_SOC <- lme(scale(OM) ~ scale(invsimpson.div), random = ~1|Replicate, data = divdat, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_SOC)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_SOC_sum <- as.data.frame(summary(lmtest_SOC)$tTable)
divdat$predlmSOC[which(is.na(divdat$OM)==F)] <- predict(lmtest_SOC, level=0)



df <- data.frame(statistic = c("intercept", "slope", "t", "df", "p"),
                 value_POXC = c(lmtest_POXC_sum$Value[1], lmtest_POXC_sum$Value[2], lmtest_POXC_sum$`t-value`[2], lmtest_POXC_sum$DF[2], lmtest_POXC_sum$`p-value`[2]),
                 value_SOC = c(lmtest_SOC_sum$Value[1], lmtest_SOC_sum$Value[2], lmtest_SOC_sum$`t-value`[2], lmtest_SOC_sum$DF[2], lmtest_SOC_sum$`p-value`[2]))


write.csv(df, "Model-output/effects_diversity/Cloutier_lmm-carbon_amf.csv")





p <- ggplot(divdat, aes(x=invsimpson.div, y=OM)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="AMF Inverse Simpson diversity", y=expression(paste("Soil organic matter (%)")), fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cover)) +
  theme_bw()  
p

q <- ggplot(divdat, aes(x=invsimpson.div, y=POXC)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="AMF Inverse Simpson diversity", y=expression(paste("Permanganate oxidizable carbon (mg kg"^-1,")")), fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cover)) +
  theme_bw()  +
  theme(axis.title=element_text(size=10))
q

  

r <- ggpubr::ggarrange(p+ ggpubr::rremove("xlab"),q+ ggpubr::rremove("xlab"),
                       nrow=1, common.legend=TRUE, legend="right")
s <- ggpubr::annotate_figure(r, bottom = ggpubr::text_grob("AMF Inverse Simpson diversity"))
ggpubr::ggexport(s, filename="Figures/regressions/Cloutier-invsimpson.div_fung_AMF.pdf", height=3.25,width=8)


