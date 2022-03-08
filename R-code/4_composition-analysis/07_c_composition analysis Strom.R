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



path1 <- "Raw-data/study-sequences/Strom-2020/seqs" # path for fungi




################  FUNGI ################  
# asv taxonomy
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqtabs/seqtab3_eukrem.csv"))


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.1", by.y = "X.1", 
                  all.x = TRUE, all.y = FALSE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.1 and X.x columns in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- 19:114
divmat <- t(taxabund[,plots])
colnames(divmat) <- taxabund$X.x
rownames(divmat) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat)
range(colSums(divmat)) # check that all asvs have at least 10 reads


# dataframe for diversity estimates and treatments
meta <- read.csv("Raw-data/study-sequences/Strom-2020/SraRunTable.csv")
meta2 <- meta[which(meta$Experiment %in% rownames(divmat)),]
meta2$Yield <- meta2$SoybeanYield 
meta2$Yield[which(meta2$Crop=="Corn")] <- meta2$CornYield[which(meta2$Crop=="Corn")]
divdat <- data.frame(Cropping.system = meta2$crop_rotation,
                     Season = meta2$season,
                     Year = meta2$year,
                     Replicate = meta2$block, 
                     ID = meta2$Experiment, 
                     Yield =  meta2$Yield,
                     TOC =  meta2$TOC)
divdat$Cropping.system0 <- rep("Corn", dim(divdat)[1])
divdat$Cropping.system0[which(divdat$Cropping.system=="Ca")] <- "(Corn)-soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="Sa")] <- "Corn-(soybean)"
divdat$Cropping.system0[which(divdat$Cropping.system=="Ss")] <- "Soybean"
divdat$Cropping.system0 <- as.factor(divdat$Cropping.system0)
divdat$Cropping.system0 <- factor(divdat$Cropping.system0, levels(divdat$Cropping.system0)[c(2,1,4,3)])
divdat$Cropping.system <- divdat$Cropping.system0
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$ID2 <- paste(divdat$Cropping.system, divdat$Replicate, sep="_")
divdat$Season <- as.factor(divdat$Season)
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(3,2,1)])
levels(divdat$Season) <- c("Spring", "Summer", "Fall")
divdat$Year <- as.factor(divdat$Year)


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
myfileENV <- myfileENV[order(myfileENV$ID),]
rownames(myfileENV) <- myfileENV$ID

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

#  data
ENV <- myfileENV
OTU <- myfileOTU


#  analysis
rownames(t(OTU))==rownames(ENV) # check


vdat <- vegdist(t(OTU), method="jaccard")
mod<-dbrda(vdat ~ TOC + Condition(Replicate), data=ENV, na.action="na.omit")

h<-how(plots=Plots(strata=ENV$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=999))
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/strom_anova.csv")
summod <- summary(mod)
write.csv(round(as.data.frame(summod[c("tot.chi", "partial.chi", "constr.chi", "unconst.chi")]),3), "Model-output/db-rda/strom_chi.csv")

bis <- as.data.frame(summod$biplot)
rownames(bis)
rownames(bis) <- c("SOC ***")
bis$color <- c("black")


#  plot

cropgroups<-ENV$Cropping.system
ssngroups <-ENV$Season
yeargroups <-ENV$Year
cols <- RColorBrewer::brewer.pal(n = 6, name = "YlGnBu")[c(3:6)]
pts <- c(21,22,23)

pdf("Figures/communities/Strom_rda.pdf", height=4, width=4)
par(mar=c(1.2,1.95,2,0.1))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=cols[cropgroups],
     pch=pts[ssngroups], 
     col=c(NA, "black")[yeargroups], cex=1.25, xlim=c(-2.3, 2.4), 
     las=1, mgp = c(2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cropping system", -2.3,-0.5, legend=levels(ENV$Cropping.system), col="white",
       pt.bg=cols, pch=21, pt.cex=1.25, cex=0.75)
legend(title="Season", 1.35,1.6, legend=levels(ENV$Season), 
       pch=pts, pt.cex=1, cex=0.75, bg="white")
legend(title="Year", 1.6, 0.8, legend=c("2015", "2016"), pt.bg = "gray50",
       col=c(NA, "black"), pch=21, pt.cex=1, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=0, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.4, 0, rownames(bis), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("C. Strom ",italic("et al.")," (2020)")), side=3, adj=0, line=0.25, cex=1, font=4)
dev.off()










#######################################





### permanova 
colnames(divmat_pa)[1:5]; rownames(divmat_pa)[1:5] # asvs as columns, samples as rows
colnames(divdat)[1:2]; rownames(divdat)[1:5] # treatment levels as columns, samples as rows
# permanova model
permmod <- adonis(divmat_pa ~ Cropping.system, 
                  data=divdat, 
                  method="bray", 
                  strata=divdat$Replicate) 
permmod2 <- round(as.data.frame(permmod$aov.tab), 2)
# append to a mod.effects dataframe
mod.effects <- data.frame(effect=c("Cropping.system")) # main effect
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$result <- ssp[1]
write.csv(mod.effects, "Model-output/Strom_permanova_fun.csv")
# examine treatment effects








