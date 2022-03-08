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
divdat <- read.csv( "Processed-data/Strom-diversity.csv")
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(2,1,4,3)])
divdat$Year <- as.factor(divdat$Year)
divdat$Season <- as.factor(divdat$Season)
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(3,2,1)])
divdat$Replicate <- as.factor(divdat$Replicate)


colorvec <- c("burlywood4", "darkolivegreen4", "burlywood4", "darkolivegreen4")



################  MICROBIAL DIVERSITY RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Year",
                                   "Season",
                                   "Cropping.system:Year",
                                   "Cropping.system:Season",
                                   "Year:Season",
                                   "Cropping.system:Year:Season")) 
options("scipen"=100, "digits"=4)




## inverse Simpson's diversity ## 
divdat$invsimpson.div <- as.numeric(divdat$invsimpson.div)
hist((divdat$invsimpson.div))
# mixed effects model
fit <- lme(invsimpson.div ~ Cropping.system*Year*Season, random = ~1|Replicate, data = divdat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(divdat[which(is.na(divdat$invsimpson.div)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Year, Season)%>% 
  summarise_at(vars(invsimpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$invsimpson.div <- ssp

# marginal interaction
test1 <- emmeans(fit, ~Cropping.system*Season)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system, testlet$Season),] # no differences

# visualize
### tukey & visualize
# compute least-squares means
test1 <- emmeans(fit, ~Cropping.system)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Cropping.system),]
p <- ggplot(divdat, aes(x=Cropping.system, y=invsimpson.div, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Fungal Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(divdat$invsimpson.div)) + c(0,11))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1,2,3,4), y = rep(125, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Strom_fun_cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Strom_diversity_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

numcomp <- 2
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "(Corn)-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cropping.system==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 2
comp <- 2
Lt <- "Corn-(soybean)" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- divdat$invsimpson.div[which(divdat$Cropping.system==Lt)] # diversity treatment mean
Yc <- divdat$invsimpson.div[which(divdat$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


write.csv(d, "Model-output/effects_diversity/Strom_diversity_RR.csv")








################  CARBON RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system",
                                   "Year",
                                   "Season",
                                   "Cropping.system:Year",
                                   "Cropping.system:Season",
                                   "Year:Season",
                                   "Cropping.system:Year:Season")) 
options("scipen"=100, "digits"=4)




## SOC ## 
divdat$SOC <- as.numeric(divdat$SOC)
hist((divdat$SOC)^3)
# mixed effects model
fit <- lme((SOC)^3 ~ Cropping.system*Year*Season, random = ~1|Replicate, data = divdat)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(divdat[which(is.na(divdat$SOC)==FALSE),],resid)
temp = resid%>%
  group_by(Cropping.system, Year, Season)%>% 
  summarise_at(vars(SOC),funs(mean=mean(., na.rm=T),Std=sd(.)))
max(temp$Std)/min(temp$Std) # check that max/min variation across groups is less than five-fold
# anova table
mod.out <- as.data.frame(car::Anova(fit)) # default SS for Anova is type-III
mod.out <- mod.out %>% mutate_if(is.numeric, round, digits=3)
# append to mod.effects dataframe
ssp <- paste0("Chisq=", mod.out$`Chisq`, ", Df=", mod.out$`Df`, ", p=", mod.out$`Pr(>Chisq)`)
mod.effects$SOC <- ssp

#  interaction
test1 <- emmeans(fit, ~Cropping.system*Year)
# compact letter display
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed=T)
testlet <- testlet[order(testlet$Year, testlet$Cropping.system),] 

p <- ggplot(divdat, aes(x=Cropping.system, y=SOC, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(divdat$SOC)) + c(0,0.1))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("burlywood4", "darkolivegreen4", "burlywood4", "darkolivegreen4")) +
  guides(fill=FALSE) +
  geom_text(x = rep(c(1,2,3,4),2), y = rep(3.3, 8), aes(label = .group), data = testlet) +
  facet_grid(.~Year)
p
ggpubr::ggexport(p, height=3, width=9, filename = "Figures/effect-plots/Strom_SOC_cropsys-yr.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Strom_carbon_anova.csv")




## effect sizes for carbon response to crop diversification ##
# use log-response ratio 
numcomp <- 4
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp),
                Year=rep(NA, numcomp))

# comparison 1
comp <- 1
Year <- 2015
Lt <- "(Corn)-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- divdat$SOC[which(divdat$Cropping.system==Lt & divdat$Year==Year)] # diversity treatment mean
Yc <- divdat$SOC[which(divdat$Cropping.system==Lc & divdat$Year==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- Year

# comparison 2
comp <- 2
Year <- 2015
Lt <- "Corn-(soybean)" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- divdat$SOC[which(divdat$Cropping.system==Lt & divdat$Year==Year)] # diversity treatment mean
Yc <- divdat$SOC[which(divdat$Cropping.system==Lc & divdat$Year==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- Year


# comparison 3
comp <- 3
Year <- 2016
Lt <- "(Corn)-soybean" # label for diversity treatment
Lc <- "Corn" # label for diversity control
Yt <- divdat$SOC[which(divdat$Cropping.system==Lt & divdat$Year==Year)] # diversity treatment mean
Yc <- divdat$SOC[which(divdat$Cropping.system==Lc & divdat$Year==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- Year

# comparison 4
comp <- 4
Year <- 2016
Lt <- "Corn-(soybean)" # label for diversity treatment
Lc <- "Soybean" # label for diversity control
Yt <- divdat$SOC[which(divdat$Cropping.system==Lt & divdat$Year==Year)] # diversity treatment mean
Yc <- divdat$SOC[which(divdat$Cropping.system==Lc & divdat$Year==Year)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d$Year[comp] <- Year



write.csv(d, "Model-output/effects_diversity/Strom_carbon_RR.csv")







################  MICROBIAL DIVERSITY ASSOCIATION WITH SOIL CARBON  ################  

# scale x and y variables to get standardized beta coefficients


lmtest_SOC <- lme(scale(SOC) ~ scale(invsimpson.div), random = ~1|Replicate, data = divdat, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_SOC)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_SOC_sum <- as.data.frame(summary(lmtest_SOC)$tTable)



df <- data.frame(statistic = c("intercept", "slope", "t", "df", "p"),
                 value_SOC = c(lmtest_SOC_sum$Value[1], lmtest_SOC_sum$Value[2], lmtest_SOC_sum$`t-value`[2], lmtest_SOC_sum$DF[2], lmtest_SOC_sum$`p-value`[2]))


write.csv(df, "Model-output/effects_diversity/Strom_lmm-carbon.csv")







p <- ggplot(divdat, aes(x=invsimpson.div, y=SOC)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="Fungal Inverse Simpson diversity", y=expression(paste("Soil organic carbon (%)")), fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cropping.system)) +
  theme_bw()  +
  scale_fill_manual(values=c("burlywood","burlywood4", "darkolivegreen2","darkolivegreen4")) +
  theme(axis.title=element_text(size=12))
p
ggpubr::ggexport(p, height=3.2, width=4.75, filename = "Figures/regressions/Strom_SOC-invsimpson.div_fun.pdf")


