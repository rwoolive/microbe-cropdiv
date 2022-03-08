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




# Study metadata
meta <- read.csv("Raw-data/study-sequences/Ashworth/Ash-Metadata.csv")
meta %>% group_by(SequenceLet, Yr, Location) %>% summarize(mean=mean(Yield))
meta$yield.crop <- rep("corn", dim(meta)[1])
meta$yield.crop[which(meta$SequenceLet=="S>S>S>S")] <- "soybean"
meta$yield.crop[which(meta$SequenceLet=="C>S>C>S" & meta$Yr==1)] <- "soybean"
meta$yield.crop[which(meta$SequenceLet=="Ct>Ct>Ct>Ct")] <- "cotton"

# Study microbial diversity data
divdat <- read.csv("Processed-data/Ashworth_both-sites-diversity.csv")
divdat_mtrec <- read.csv("Processed-data/Ashworth_MTREC-diversity.csv")
divdat_recm <- read.csv("Processed-data/Ashworth_recm-diversity.csv")


# combine data into one dataset
dat <- merge(meta, divdat, by.x="Sequence", by.y="ID")
dat_mtrec <- merge(meta, divdat_mtrec, by.x="Sequence", by.y="ID")
dat_recm <- merge(meta, divdat_recm, by.x="Sequence", by.y="ID")

dat$Cropping.system <- as.factor(dat$Cropping.system)
dat$Cropping.system <- factor(dat$Cropping.system, levels(dat$Cropping.system)[c(1,3,2)])
dat$Cover <- as.factor(dat$Cover.x)
dat$Replicate <- as.factor(dat$Replicate)
dat$Yr <- as.factor(dat$Yr.x)

dat_mtrec$Cropping.system <- as.factor(dat_mtrec$Cropping.system)
dat_mtrec$Cropping.system <- factor(dat_mtrec$Cropping.system, levels(dat_mtrec$Cropping.system)[c(1,3,2)])
dat_mtrec$Cover <- as.factor(dat_mtrec$Cover.x)
dat_mtrec$Replicate <- as.factor(dat_mtrec$Replicate)
dat_mtrec$Yr <- as.factor(dat_mtrec$Yr.x)

dat_recm$Cropping.system <- as.factor(dat_recm$Cropping.system)
dat_recm$Cropping.system <- factor(dat_recm$Cropping.system, levels(dat_recm$Cropping.system)[c(1,4,3,2)])
dat_recm$Cover <- as.factor(dat_recm$Cover.x)
dat_recm$Replicate <- as.factor(dat_recm$Replicate)
dat_recm$Yr <- as.factor(dat_recm$Yr.x)



# remove cotton monoculture and chicken litter treatments
dat <- dat[-which(dat$Cover=="Litter"),]
dat$Cover <- factor(dat$Cover, levels(dat$Cover)[c(1,3,4)])

dat_mtrec <- dat_mtrec[-which(dat_mtrec$Cover=="Litter"),]
dat_mtrec$Cover <- factor(dat_mtrec$Cover, levels(dat_mtrec$Cover)[c(1,3,4)])


dat_recm <- dat_recm[-which(dat_recm$Cover=="Litter" | dat_recm$Cropping.system=="Cotton"),]
dat_recm$Cover <- factor(dat_recm$Cover, levels(dat_recm$Cover)[c(1,3,4)])
dat_recm$Cropping.system <- factor(dat_recm$Cropping.system, levels(dat_recm$Cropping.system)[c(1,2,4)])



colorvec <- c("burlywood4", rep("darkolivegreen3", 2))
colorvec2 <- c(rep("burlywood4", 2), "darkolivegreen3")



# TEST FOR INTERACTIONS BETWEEN DIVERSIFICATION AND OTHER FACTORS
dat$invsimpson.div <- as.numeric(dat$invsimpson.div)
dat$Location <- as.factor(dat$Location.x)
hist((dat$invsimpson.div))
# mixed effects model
fit <- lme(invsimpson.div ~ Location*Cropping.system*Cover*Yr, random = ~1|Replicate/Yr/Location/Cropping.system, data = dat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat[which(is.na(dat$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Location,Cropping.system, Cover, Yr)%>% 
  summarise_at(vars(invsimpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)


################  MICROBIAL DIVERSITY RESPONSE TO CROP DIVERSIFICATION  ################  

### MTREC ### 

### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Cover",
                                   "Yr",
                                   "Cropping.system:Cover",
                                   "Cropping.system:Yr",
                                   "Cover:Yr",
                                   "Cropping.system:Cover:Yr")) 
options("scipen"=100, "digits"=4)




## inverse Simpson's diversity ## 
dat_mtrec$invsimpson.div <- as.numeric(dat_mtrec$invsimpson.div)
hist((dat_mtrec$invsimpson.div))
# mixed effects model
fit <- lme(invsimpson.div ~ Cropping.system*Cover*Yr, random = ~1|Replicate/Cropping.system, data = dat_mtrec)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat_mtrec[which(is.na(dat_mtrec$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Cover, Yr)%>% 
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
p <- ggplot(dat_mtrec, aes(x=Cover, y=invsimpson.div, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cover", y=expression(paste("Bacterial Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat_mtrec$invsimpson.div)) + c(0,100))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(1200, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_mtrec-cover.pdf")

test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(dat_mtrec, aes(x=Cropping.system, y=invsimpson.div, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Bacterial Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat_mtrec$invsimpson.div)) + c(0,100))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec2) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(1200, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_mtrec-cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Ashworth_diversity-mtrec_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

## MTREC, cover comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Vetch" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cover==Lt)] # diversity treatment values
Yc <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cover==Lc)] # diversity control values

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Wheat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cover==Lt)] # diversity treatment values
Yc <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cover==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Ashworth_diversity-mtrec-cover_RR.csv")


## MTREC, cropsys comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cropping.system==Lc)] # diversity control values

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_mtrec$invsimpson.div[which(dat_mtrec$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



write.csv(d, "Model-output/effects_diversity/Ashworth_diversity-mtrec-cropsys_RR.csv")








### recm ### 

### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Cover",
                                   "Yr",
                                   "Cropping.system:Cover",
                                   "Cropping.system:Yr",
                                   "Cover:Yr",
                                   "Cropping.system:Cover:Yr")) 
options("scipen"=100, "digits"=4)




## inverse Simpson's diversity ## 
dat_recm$invsimpson.div <- as.numeric(dat_recm$invsimpson.div)
hist((dat_recm$invsimpson.div))
# mixed effects model
fit <- lme(invsimpson.div ~ Cropping.system*Cover*Yr, random = ~1|Replicate/Cropping.system/Cover, data = dat_recm)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat_recm[which(is.na(dat_recm$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Cover, Yr)%>% 
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
test1 <- emmeans(fit, ~Cover, by="Yr")
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet_yr <- testlet[order(testlet$Yr, testlet$Cover),]
p <- ggplot(dat_recm, aes(x=Cover, y=invsimpson.div, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cover", y=expression(paste("Bacterial Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat_recm$invsimpson.div)) + c(0,200))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  facet_grid(Yr~.) +
  guides(fill=FALSE) +
  geom_text(x = rep(c(1:3), 2), y = rep(1100, 3*2), aes(label = .group), data = testlet_yr)
p
ggpubr::ggexport(p, height=6, width=4, filename = "Figures/effect-plots/Ashworth_recm-cover-yr.pdf")

test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(dat_recm, aes(x=Cropping.system, y=invsimpson.div, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Bacterial Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat_recm$invsimpson.div)) + c(0,200))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec2) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(1100, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_recm-cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Ashworth_diversity-recm_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

## recm, cover comparisons
numcomp <- 4

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp),
                Year=rep(NA,numcomp))

# comparison 1
comp <- 1
Lt <- "Vetch" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Year <- 1
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lt & dat_recm$Yr==Year)] # diversity treatment mean
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lc & dat_recm$Yr==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- paste(Year)


# comparison 2
comp <- 2
Lt <- "Vetch" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Year <- 2
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lt & dat_recm$Yr==Year)] # diversity treatment mean
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lc & dat_recm$Yr==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- paste(Year)



# comparison 3
comp <- 3
Lt <- "Wheat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Year <- 1
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lt & dat_recm$Yr==Year)] # diversity treatment mean
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lc & dat_recm$Yr==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- paste(Year)


# comparison 4
comp <- 4
Lt <- "Wheat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Year <- 2
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lt & dat_recm$Yr==Year)] # diversity treatment mean
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cover==Lc & dat_recm$Yr==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- paste(Year)


write.csv(d, "Model-output/effects_diversity/Ashworth_diversity-recm-cover-yr_RR.csv")


## recm, cropsys comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- dat_recm$invsimpson.div[which(dat_recm$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_recm$invsimpson.div[which(dat_recm$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)




write.csv(d, "Model-output/effects_diversity/Ashworth_diversity-recm-cropsys_RR.csv")








################  CARBON RESPONSE TO CROP DIVERSIFICATION  ################  

### MTREC ### 


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Cover",
                                   "Yr",
                                   "Cropping.system:Cover",
                                   "Cropping.system:Yr",
                                   "Cover:Yr",
                                   "Cropping.system:Cover:Yr")) 
options("scipen"=100, "digits"=4)





## C.s ## 
dat_mtrec$C.s <- as.numeric(dat_mtrec$C.s)
hist((dat_mtrec$C.s))
# mixed effects model
fit <- lme((C.s) ~ Cropping.system*Cover*Yr, random = ~1|Replicate/Cropping.system, data = dat_mtrec, na.action="na.omit")
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat_mtrec[which(is.na(dat_mtrec$C.s)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Cover, Yr)%>% 
  summarise_at(vars(C.s),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$C.s <- ssp

# interaction but no effects of diversity
test1 <- emmeans(fit, ~Cropping.system*Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system, testlet$Cover),] # no differences
testlet <- testlet[order(testlet$Cover, testlet$Cropping.system),] # no differences

# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cover),]
p <- ggplot(dat_mtrec, aes(x=Cover, y=C.s, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cover", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(dat_mtrec$C.s)) + c(0,0.13))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(2.2, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_mtrec-carbon_cover.pdf")

test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(dat_mtrec, aes(x=Cropping.system, y=C.s, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(dat_mtrec$C.s)) + c(0,0.13))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec2) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(2.2, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_mtrec-carbon_cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Ashworth_carbon-mtrec_anova.csv")




## effect sizes for carbon response to crop diversification ##
# use log-response ratio 

## MTREC, cover comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Vetch" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_mtrec$C.s[which(dat_mtrec$Cover==Lt)] # diversity treatment values
Yc <- dat_mtrec$C.s[which(dat_mtrec$Cover==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Wheat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_mtrec$C.s[which(dat_mtrec$Cover==Lt)] # diversity treatment values
Yc <- dat_mtrec$C.s[which(dat_mtrec$Cover==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Ashworth_carbon-mtrec-cover_RR.csv")


## MTREC, cropsys comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- dat_mtrec$C.s[which(dat_mtrec$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_mtrec$C.s[which(dat_mtrec$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- dat_mtrec$C.s[which(dat_mtrec$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_mtrec$C.s[which(dat_mtrec$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Ashworth_carbon-mtrec-cropsys_RR.csv")





### recm ### 


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Cover",
                                   "Yr",
                                   "Cropping.system:Cover",
                                   "Cropping.system:Yr",
                                   "Cover:Yr",
                                   "Cropping.system:Cover:Yr")) 
options("scipen"=100, "digits"=4)





## C.s ## 
dat_recm$C.s <- as.numeric(dat_recm$C.s)
hist((dat_recm$C.s))
# mixed effects model
fit <- lme((C.s) ~ Cropping.system*Cover*Yr, random = ~1|Replicate/Cropping.system, data = dat_recm, na.action="na.omit")
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat_recm[which(is.na(dat_recm$C.s)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Cover, Yr)%>% 
  summarise_at(vars(C.s),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$C.s <- ssp


# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cover)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cover),]
p <- ggplot(dat_recm, aes(x=Cover, y=C.s, fill=Cover)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cover", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(dat_recm$C.s)) + c(0,0.13))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(2.6, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_recm-carbon_cover.pdf")

test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(dat_recm, aes(x=Cropping.system, y=C.s, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(dat_recm$C.s)) + c(0,0.13))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec2) +
  guides(fill=FALSE) +
  geom_text(x = c(1:3), y = rep(2.6, 3), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=4, filename = "Figures/effect-plots/Ashworth_recm-carbon_cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Ashworth_carbon-recm_anova.csv")






## effect sizes for carbon response to crop diversification ##
# use log-response ratio 

## recm, cover comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Vetch" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_recm$C.s[which(dat_recm$Cover==Lt)] # diversity treatment values
Yc <- dat_recm$C.s[which(dat_recm$Cover==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Wheat" # label for diversity treatment
Lc <- "Fallow" # label for diversity control
Yt <- dat_recm$C.s[which(dat_recm$Cover==Lt)] # diversity treatment values
Yc <- dat_recm$C.s[which(dat_recm$Cover==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Ashworth_carbon-recm-cover_RR.csv")


## recm, cropsys comparisons
numcomp <- 2

d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- dat_recm$C.s[which(dat_recm$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_recm$C.s[which(dat_recm$Cropping.system==Lc)] # diversity control values


d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- dat_recm$C.s[which(dat_recm$Cropping.system==Lt)] # diversity treatment values
Yc <- dat_recm$C.s[which(dat_recm$Cropping.system==Lc)] # diversity control values

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Ashworth_carbon-recm-cropsys_RR.csv")










################  MICROBIAL DIVERSITY ASSOCIATION WITH SOIL CARBON  ################  
# scale x and y variables to get standardized beta coefficients

lmtest_mtrec <- lme(scale(C.s) ~ scale(invsimpson.div), random = ~1|Replicate, data = dat_mtrec, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_mtrec)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_mtrec_sum <- as.data.frame(summary(lmtest_mtrec)$tTable)
dat_mtrec$predlm[which(is.na(dat_mtrec$C.s)==F)] <- predict(lmtest_mtrec, level=0)


lmtest_recm <- lme(scale(C.s) ~ scale(invsimpson.div), random = ~1|Replicate, data = dat_recm)
# check redisuals
resid <- residuals(lmtest_recm)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_recm_sum <- as.data.frame(summary(lmtest_recm)$tTable)
dat_recm$predlm[which(is.na(dat_recm$C.s)==F)] <- predict(lmtest_recm, level=0)


df <- data.frame(statistic = c("intercept", "slope", "t", "df", "p"),
           value_mtrec = c(lmtest_mtrec_sum$Value[1], lmtest_mtrec_sum$Value[2], lmtest_mtrec_sum$`t-value`[2], lmtest_mtrec_sum$DF[2], lmtest_mtrec_sum$`p-value`[2]),
           value_recm = c(lmtest_recm_sum$Value[1], lmtest_recm_sum$Value[2], lmtest_recm_sum$`t-value`[2], lmtest_recm_sum$DF[2], lmtest_recm_sum$`p-value`[2]))

write.csv(df, "Model-output/effects_diversity/Ashworth_lmm-carbon.csv")





# plots
p <- ggplot(dat_mtrec, aes(x=scale(invsimpson.div), y=scale(C.s))) +
  geom_line(aes(y=predlm), size=1, color="black", linetype="dashed") +
  labs(x="Bacterial Inverse Simpson diversity", y=expression(paste("Soil organic carbon (%)")), fill="Cropping system", subtitle="A. Ashworth et al. (2017) @ MTREC") +
  geom_point(aes(fill=Cropping.system, shape=Cover, size=Yr)) +
  scale_fill_manual(values=c("burlywood4", "burlywood", "darkolivegreen4")) +
  scale_shape_manual(values=c(22,23,24)) +
  scale_size_manual(values=c(2,3)) +
  theme_bw()  + 
  guides(fill = guide_legend(override.aes=list(shape=21))) 
p
ggpubr::ggexport(p, height=3.5, width=6, filename = "Figures/regressions/Ashworth_soc-invsimpson.div_bac_MTREC.pdf")



# plots
p <- ggplot(dat_recm, aes(x=invsimpson.div, y=C.s)) +
  geom_line(aes(y=predlm), size=1, color="black", linetype="dashed") +
  labs(x="Bacterial Inverse Simpson diversity", y=expression(paste("Soil organic carbon (%)")), fill="Cropping system", subtitle="B. Ashworth et al. (2017) @ RECM") +
  geom_point(aes(fill=Cropping.system, shape=Cover, size=Yr)) +
  scale_fill_manual(values=c("burlywood4", "burlywood", "darkolivegreen4")) +
  scale_shape_manual(values=c(22,23,24)) +
  scale_size_manual(values=c(2,3)) +
  theme_bw()  + 
  guides(fill = guide_legend(override.aes=list(shape=21))) 
p
ggpubr::ggexport(p, height=3.5, width=6, filename = "Figures/regressions/Ashworth_soc-invsimpson.div_bac_RECM.pdf")



