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







################  Bacteria ################  
path1 <- "Raw-data/study-sequences/Ashworth/All" # path for bacteria

authorM <- "Ashworth"
sampdateM <- c("year-1", "year-2") # sampling dates
locationM <- c("mtrec","recm") # sampling locations
crop.seqM <- c("cccc","cscs","ssss","ctctctct") # crop sequences
regionM <- "16S"




# asv taxonomy
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqtab/3seqtab.csv"))


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.1", by.y = "X.1", 
                  all.x = TRUE, all.y = TRUE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.1 and X.x columns in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- which(substr(colnames(taxabund), 2, 4) == regionM) # columns that contain abundances
divmat <- taxabund[,plots]
divmat2 <- t(divmat)
colnames(divmat2) <- taxabund$X.x
rownames(divmat2) <- colnames(taxabund)[plots]
rownames(divmat2) <- gsub("X", "", rownames(divmat2))
divmat <- as.data.frame(divmat2)
range(colSums(divmat)) # check that all asvs have at least 10 reads


# dataframe for diversity estimates and treatments
meta <- read.csv("Raw-data/study-sequences/Ashworth/Ash-Metadata.csv")
meta$Sequence2 <- paste0(regionM, ".", meta$Location, ".", meta$SequenceLet, ".", meta$Cover, 
                         ".Rep", meta$Rep, ".Yr", meta$Yr)
meta$Sequence2 <- gsub(">", ".", meta$Sequence2)
meta2 <- meta[which(meta$Sequence2 %in% rownames(divmat)),]
write.csv(meta2, "Raw-data/study-sequences/Ashworth/Ash-Metadata-analyzed.csv")
rownames(meta2) <- meta2$Sequence2
meta2 <- meta2[rownames(divmat),]
divdat <- data.frame(Location = as.factor(meta2$Location),
                     Cropping.system0 = as.factor(meta2$SequenceLet),
                     Cover = as.factor(meta2$Cover),
                     Replicate = as.factor(meta2$Rep), 
                     Yr = as.factor(meta2$Yr),
                     ID = meta2$Sequence, 
                     SOC = meta2$C.s,
                     Yield = meta2$Yield)
divdat$Cropping.system <- revalue(divdat$Cropping.system0, c("C>C>C>C"="Corn", "C>S>C>S"="Corn-soybean", "Ct>Ct>Ct>Ct"="Cotton", "S>S>S>S"="Soybean"))
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])
divdat$ID2 <- paste0("16S.", meta2$Location, ".", meta2$SequenceLet, ".", meta2$Cover, ".Rep", meta2$Rep, ".Yr", meta2$Yr)
divdat$ID2 <- gsub(">", ".", divdat$ID2)

table(divdat$Cropping.system, divdat$Cover, divdat$Location, divdat$Yr)

rownames(divmat)==divdat$ID2 # check that both dataframes are ordered by ID









# divmat_pa: matrix with ASVs as columns and plots as rows (proportional abundance)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_pa <- proportions(as.matrix(divmat),1)
range(rowSums(divmat_pa)) # check to see that all rows sum up to 1
divmat_pa <- data.frame(divmat_pa)




# look at summary of reads for each sample
summary(rowSums(divmat))
summary(rowSums(divmat_pa))





######
## RDA ordination 

# proportional abundance OTU file with ASVs as rows and samples as columns
myfileOTU <- divmat_pa
head(myfileOTU)[1:5,1:5]
myfileOTU <- t(myfileOTU)
head(myfileOTU)[1:5,1:5]
myfileOTU <- myfileOTU[,order(colnames(myfileOTU))]

# Environment file with samples as rows and properties as columns
myfileENV <- divdat
myfileENV <- myfileENV[order(myfileENV$ID2),]
rownames(myfileENV) <- myfileENV$ID2

# make sure that both datasets are ordered in the same way
check <- data.frame(myfileOTU = colnames(myfileOTU), 
           myfileENV = rownames(myfileENV))
check$myfileOTU==check$myfileENV

library(vegan)
library(ggplot2); packageVersion("ggplot2") #3.0.0
library(gdata)
library(ecodist)
library(lawstat)
library(lmerTest)
library(car)
library(multcomp)
library(PairedData)


str(myfileOTU)
str(myfileENV)

# mtrec data
ENV_mtrec <- myfileENV[which(myfileENV$Location=="MTREC"),]
ENV_mtrec$Cropping.system <- factor(ENV_mtrec$Cropping.system, levels(ENV_mtrec$Cropping.system)[c(1,2,4)])
levels(ENV_mtrec$Yr) <- c("2013", "2014")
OTU_mtrec <- myfileOTU[, which(substr(colnames(myfileOTU), 5, 8) == "MTRE")]

# recm data
ENV_recm <- myfileENV[which(myfileENV$Location=="RECM"),]
levels(ENV_recm$Yr) <- c("2013", "2014")
OTU_recm <- myfileOTU[, which(substr(colnames(myfileOTU), 5, 8) == "RECM")]



# ashworth mtrec distance-based RDA analysis
rownames(t(OTU_mtrec))==rownames(ENV_mtrec) # check

# model including treatments and SOC
vdat <- vegdist(t(OTU_mtrec), method="jaccard")
mod<-dbrda(vdat ~ SOC + Condition(Replicate), 
           data=ENV_mtrec, na.action="na.omit")

h<-how(plots=Plots(strata=ENV_mtrec$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=999)) # marginal effects are those which assume that the explanatory variable is the only one in the model
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/ashworth_mtrec_anova.csv")
summod <- summary(mod)
write.csv(round(as.data.frame(summod[c("tot.chi", "partial.chi", "constr.chi", "unconst.chi")]),3), "Model-output/db-rda/ashworth_mtrec_chi.csv")

bis <- as.data.frame(summod$biplot)
bis$color <- c("black")

# mtrec plot

cropgroups<-ENV_mtrec$Cropping.system
covgroups<-ENV_mtrec$Cover
yrgroups <-ENV_mtrec$Yr

pdf("Figures/communities/Ashworth_rda-mtrec.pdf", height=4, width=4)
par(mar=c(1.2,2.025,2,0.1))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=c("goldenrod", "green", "darkolivegreen")[cropgroups],
     pch=c(21,22,24,25)[covgroups],
     col=c(NA, "black")[yrgroups], cex=1.25, xlim=c(-3.7, 3), 
     las=1, mgp = c(2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cropping system", -3.7,2, legend=c("Corn", "Soybean", "Corn-soybean"), col="white",
       pt.bg=c("goldenrod", "green", "darkolivegreen"), pch=21, pt.cex=1.25, cex=0.75)
legend(title="Cover", -3.7,1, legend=c("Fallow", "Litter", "Vetch", "Wheat"), 
       pch=c(21,22,24,25), pt.cex=1, cex=0.75, bg = "white")
legend(title="Year", -3.7,-0.2, legend=c("2013", "2014"), pt.bg = "gray50",
       col=c(NA, "black"), pch=21, pt.cex=1, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=bis$MDS1, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.5, bis$MDS1+sign(bis$MDS1)*0.5, 
     paste(rownames(bis),"***"), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("A. Ashworth ",italic("et al.")," (2017) at MTREC")), side=3, adj=0, line=0.25, cex=1, font=4)
dev.off()




# recm analysis
rownames(t(OTU_recm))==rownames(ENV_recm) # check

vdat <- vegdist(t(OTU_recm), method="jaccard")
mod<-dbrda(vdat ~ SOC + Condition(Replicate), data=ENV_recm)

h<-how(plots=Plots(strata=ENV_recm$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=999))
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/ashworth_recm_anova.csv")
summod <- summary(mod)
bis <- as.data.frame(summod$biplot)
rownames(bis)
bis$color <- c("black")


# recm plot

cropgroups<-ENV_recm$Cropping.system
covgroups<-ENV_recm$Cover
yrgroups <-ENV_recm$Yr

pdf("Figures/communities/Ashworth_rda-recm.pdf", height=4, width=4)
par(mar=c(1.2,2.025,2,0.1))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=c("goldenrod", "green", "gray80", "darkolivegreen")[cropgroups],
     pch=c(21,22,24,25)[covgroups],
     col=c(NA, "black")[yrgroups], cex=1.25, xlim=c(-4, 2), 
     las=1, mgp = c(2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cropping system", -4,1.5, legend=c("Corn", "Soybean", "Cotton", "Corn-soybean"), col="white",
       pt.bg=c("goldenrod", "green", "gray90", "darkolivegreen"), pch=21, pt.cex=1.25, cex=0.75)
legend(title="Cover", -4,0.5, legend=c("Fallow", "Litter", "Vetch", "Wheat"), 
       pch=c(21,22,24,25), pt.cex=1, cex=0.75, bg = "white")
legend(title="Year", -4,-0.5, legend=c("2013", "2014"), pt.bg = "gray50",
       col=c(NA, "black"), pch=21, pt.cex=1, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=bis$MDS1, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.4, 0, 
     paste(rownames(bis), "*"), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("B. Ashworth ",italic("et al.")," (2017) at RECM")), side=3, adj=0, line=0.25, cex=1, font=4)

dev.off()















#######################################


### permanova 
colnames(divmat_pa)[1:5]; rownames(divmat_pa)[1:5] # asvs as columns, samples as rows
colnames(divdat)[1:2]; rownames(divdat)[1:5] # treatment levels as columns, samples as rows

divmat_pa_RECM <- divmat_pa[which(substr(rownames(divmat_pa), 5,8)=="RECM"),]
divmat_pa_MTREC <- divmat_pa[which(substr(rownames(divmat_pa), 5,8)=="MTRE"),]

mod.effects <- data.frame(effect=c("Cropping.system", "Cover", "Year", 
                                   "Cropping.system:Cover", "Cropping.system:Year",  "Cover:Year",  
                                   "Cropping.system:Cover:Year", "Residual", "Total")) #  effects


# permanova model for RECM
permmod <- adonis(divmat_pa_RECM ~ Cropping.system*Cover*Yr, 
                  data=divdat[which(divdat$Location=="RECM"),], 
                  method="bray", 
                  strata=divdat$Replicate[which(divdat$Location=="RECM")]) 
permmod2 <- round(as.data.frame(permmod$aov.tab), 2)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$result_RECM <- ssp


# permanova model for MTREC
permmod <- adonis(divmat_pa_MTREC ~ Cropping.system*Cover*Yr, 
                  data=divdat[which(divdat$Location=="MTREC"),], 
                  method="bray", 
                  permutations=divdat$Replicate[which(divdat$Location=="MTREC"),]) 
permmod2 <- round(as.data.frame(permmod$aov.tab), 2)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$result_MTREC <- ssp



write.csv(mod.effects, "Model-output/permanova/Ashworth_permanova_bac_each-site.csv")
# examine for treatment effects





