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
meta <- read.csv("Raw-data/study-sequences/Song-2018/SraRunTable_2.csv")
meta$Cropping.system.rep <- sapply(strsplit(meta$Library.Name, "-"), `[`, 1)

# Study microbial diversity data
divdat_fun <- read.csv("Processed-data/Song-diversity_fun_all-saprotrophs.csv")
divdat_fun$cs <- rep("SC", dim(divdat_fun)[1])
divdat_fun$cs[which(divdat_fun$Cropping.system=="Fallow soybean")] <- "FS"
divdat_fun$cs[which(divdat_fun$Cropping.system=="Corn-soybean")] <- "CS"
divdat_fun$cs[which(divdat_fun$Cropping.system=="Wheat-soybean")] <- "WS"
divdat_fun$Cropping.system.rep <- paste0(divdat_fun$cs, divdat_fun$Replicate)


# combine data into one dataset
dat_fun <- merge(divdat_fun, meta, by="Cropping.system.rep",
                 all.y=FALSE)

dat_fun$Cropping.system <- as.factor(dat_fun$Cropping.system)
dat_fun$Cropping.system <- factor(dat_fun$Cropping.system, levels(dat_fun$Cropping.system)[c(3,1,2,4)])
dat_fun$Replicate <- as.factor(dat_fun$Replicate)





colorvec <- c("burlywood4", "burlywood4", "darkolivegreen4", "darkolivegreen4")



################  CARBON RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system")) 
options("scipen"=100, "digits"=4)




## SOM ## 

testlet <- aggregate(SOM~Cropping.system, dat=dat_fun, FUN="mean")
testlet <- testlet[order(testlet$Cropping.system),] # no differences

p <- ggplot(dat_fun, aes(x=Cropping.system, y=SOM, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("Soil organic carbon (%)"))) +
  lims(y=c(range(na.omit(dat_fun$SOM)) + c(0,0.1))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) #+ geom_text(x = c(1,2,3,4), y = rep(3.3, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Song_SOM_cropsys.pdf")







## effect sizes for carbon response to crop diversification ##
# use log-response ratio 

numcomp <- 4
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))

# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Fallow soybean" # label for diversity control
Yt <- dat_fun$SOM[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$SOM[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)



# comparison 2
comp <- 2
Lt <- "Wheat-soybean" # label for diversity treatment
Lc <- "Fallow soybean" # label for diversity control
Yt <- dat_fun$SOM[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$SOM[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)


# comparison 3
comp <- 3
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Continuous soybean" # label for diversity control
Yt <- dat_fun$SOM[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$SOM[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)


# comparison 4
comp <- 4
Lt <- "Wheat-soybean" # label for diversity treatment
Lc <- "Continuous soybean" # label for diversity control
Yt <- dat_fun$SOM[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$SOM[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)
d[comp,4:6] <- rep(NA, 3)


write.csv(d, "Model-output/effects_diversity/Song_carbon_RR.csv")






###### FUNGI ###### 

################  MICROBIAL DIVERSITY RESPONSE TO CROP DIVERSIFICATION  ################  


### ANOVA & plotting
# Main effects are Cropping.system, year, season
# Random effect is Replicate

## Create a dataframe to store model output:
mod.effects <- data.frame(effect=c("Cropping.system")) 
options("scipen"=100, "digits"=4)




## inverse Simpson's diversity ## 
dat_fun$invsimpson.div <- as.numeric(dat_fun$invsimpson.div)
hist(sqrt(dat_fun$invsimpson.div))
# mixed effects model
fit <- lme(sqrt(invsimpson.div) ~ Cropping.system, random = ~1|Replicate, data = dat_fun)
# check redisuals
resid <- residuals(fit)
qqnorm(resid); qqline(resid); shapiro.test(resid)
# check variances
resid = cbind(dat_fun[which(is.na(dat_fun$invsimpson.div)==FALSE),],resid)
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
p <- ggplot(dat_fun, aes(x=Cropping.system, y=invsimpson.div, fill=Cropping.system)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) +
  labs(x="Cropping system", y=expression(paste("saprotroph Inverse Simpson diversity"))) +
  lims(y=c(range(na.omit(dat_fun$invsimpson.div)) + c(0,5))) +
  geom_jitter(alpha=0.75, width=0.2, height=0, shape=19, size=3) +
  theme_bw() +
  scale_fill_manual(values=colorvec) +
  guides(fill=FALSE) +
  geom_text(x = c(1,2,3,4), y = rep(57, 4), aes(label = .group), data = testlet)
p
ggpubr::ggexport(p, height=3, width=5, filename = "Figures/effect-plots/Song_fun_saprotroph_cropsys.pdf")



write.csv(t(mod.effects), "Model-output/effects_diversity/Song_fun_saprotroph-diversity_anova.csv")




## effect sizes for diversity response to crop diversification ##
# use log-response ratio 

numcomp <- 4
d <- data.frame(Comparison=rep(NA, numcomp),
                ES=rep(NA, numcomp),
                Est=rep(NA, numcomp),
                SE=rep(NA, numcomp),
                CI_lower=rep(NA, numcomp), 
                CI_upper=rep(NA, numcomp))



# comparison 1
comp <- 1
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Fallow soybean" # label for diversity control
Yt <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)



# comparison 2
comp <- 2
Lt <- "Wheat-soybean" # label for diversity treatment
Lc <- "Fallow soybean" # label for diversity control
Yt <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)

# comparison 3
comp <- 3
Lt <- "Corn-soybean" # label for diversity treatment
Lc <- "Continuous soybean" # label for diversity control
Yt <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)


# comparison 4
comp <- 4
Lt <- "Wheat-soybean" # label for diversity treatment
Lc <- "Continuous soybean" # label for diversity control
Yt <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lt)] # diversity treatment mean
Yc <- dat_fun$invsimpson.div[which(dat_fun$Cropping.system==Lc)] # diversity control mean

d$Comparison[comp] <- paste(Lt, Lc, sep=" vs. ")
d[comp,2:6] <- LRRi(Yc, Yt, confidence = 0.95)





write.csv(d, "Model-output/effects_diversity/Song_fun_saprotroph-diversity_RR.csv")











################  MICROBIAL DIVERSITY ASSOCIATION WITH SOIL CARBON  ################  

# scale x and y variables to get standardized beta coefficients


lmtest_SOM <- lme(scale(SOM) ~ scale(invsimpson.div), random = ~1|Replicate, data = dat_fun, na.action="na.omit")
# check redisuals
resid <- residuals(lmtest_SOM)
qqnorm(resid); qqline(resid); shapiro.test(resid)
lmtest_SOM_sum <- as.data.frame(summary(lmtest_SOM)$tTable)



df <- data.frame(statistic = c("intercept", "slope", "t", "df", "p"),
                 value_SOM = c(lmtest_SOM_sum$Value[1], lmtest_SOM_sum$Value[2], lmtest_SOM_sum$`t-value`[2], lmtest_SOM_sum$DF[2], lmtest_SOM_sum$`p-value`[2]))


write.csv(df, "Model-output/effects_diversity/Song_lmm-carbon_fun_saprotroph.csv")








# plots
r <- ggplot(dat_fun, aes(x=invsimpson.div, y=SOM)) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, color="black", linetype="dashed") +
  labs(x="saprotroph Inverse Simpson diversity", y=expression(paste("Soil organic matter (g kg"^-1, ")")),
       fill="Cropping system") +
  geom_point(size = 3, shape=21, aes(fill=Cropping.system)) +
  scale_fill_manual(values=c("burlywood4", "burlywood", "darkolivegreen4", "darkolivegreen2")) +
  theme_bw() 
r









s <- ggpubr::ggarrange(r, #+ ggpubr::rremove("xlab")
                       nrow=1, common.legend=TRUE, legend="right")
#s <- ggpubr::annotate_figure(r, bottom = ggpubr::text_grob("Bacterial Inverse Simpson diversity"))
ggpubr::ggexport(s, filename="Figures/regressions/Song_som-invsimpson.div_fun_saprotroph.pdf", height=3.25,width=5)





