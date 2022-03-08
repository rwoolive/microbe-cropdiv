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



path1 <- "Raw-data/study-sequences/Song-2018/seqs/bac" # path for bacteria
path2 <- "Raw-data/study-sequences/Song-2018/seqs/fun" # path for fungi




################  Bacteria ################  
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
plots <- 17:28
divmat <- t(taxabund[,plots])
colnames(divmat) <- taxabund$X.x
rownames(divmat) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat)
range(colSums(divmat)) # check that all asvs have at least 10 reads

# dataframe for diversity estimates and treatments
divdat <- data.frame(Cropping.system = substr(as.factor(rownames(divmat)), 1, 2),
                     Replicate = substr(as.factor(rownames(divmat)), 3, 3))
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(2,3,1,4)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Cropping.system.rep <- paste0(divdat$Cropping.system, divdat$Replicate)

# metadata
meta <- read.csv("Raw-data/study-sequences/Song-2018/SraRunTable_2.csv")
meta$Cropping.system.rep <- sapply(strsplit(meta$Library.Name, "-"), `[`, 1)
meta <- meta[,c("Cropping.system.rep", "SOM")]

divdat <- merge(divdat, meta, by="Cropping.system.rep")
divdat$names <- paste0(divdat$Cropping.system.rep,".bacterial")

divdat$names == rownames(divmat)
divdat$ID <- paste0(divdat$Cropping.system.rep, ".bacterial")



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

ENV$SOM2 <- ENV$SOM + rnorm(12,0.0001,0.0001)
vdat <- vegdist(t(OTU), method="jaccard")
mod<-dbrda(vdat ~ SOM +Condition(Replicate), data=ENV, na.action="na.omit")

h<-how(plots=Plots(strata=ENV$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=9999))
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/song_anova_bac.csv")
summod <- summary(mod)
write.csv(round(as.data.frame(summod[c("tot.chi", "partial.chi", "constr.chi", "unconst.chi")]),3), "Model-output/db-rda/song_chi_bac.csv")

bis <- as.data.frame(summod$biplot)
rownames(bis)
rownames(bis) <- c("SOC **")
bis$color <- c("black")


#  plot

cropgroups<-ENV$Cropping.system
levels(ENV$Cropping.system) <- c("Fallow soybean", "Continuous soybean", "Corn-soybean", "Wheat-soybean")
cols <- RColorBrewer::brewer.pal(n = 4, name = "Paired")

pdf("Figures/communities/Song_rda_bac.pdf", height=4, width=4)
par(mar=c(1.2,1.95,2,0.3))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=cols[cropgroups],
     pch=21, cex=1.25, xlim=c(-1, 1.55), 
     las=1, mgp = c(2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cropping system", 0.5,1.3, legend=levels(ENV$Cropping.system), col="black",
       pt.bg=cols, pch=21, pt.cex=1.25, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=0, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.25, 0, rownames(bis), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("D. Song ",italic("et al.")," (2018)")), side=3, adj=0, line=0.25, cex=1, font=4)

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
write.csv(mod.effects, "Model-output/Song-2018_permanova_bac.csv")
# examine for treatment effects














################  Fungi ################  
# asv taxonomy
taxa <- read.csv(paste0(path2, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path2, "/seqtabs/seqtab3_eukrem.csv"))


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.5", by.y = "X.1", 
                  all.x = TRUE, all.y = TRUE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.1 and X.x columns in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- 23:34
divmat <- t(taxabund[,plots])
colnames(divmat) <- taxabund$X.x
rownames(divmat) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat)
range(colSums(divmat)) # check that all asvs have at least 10 reads

# dataframe for diversity estimates and treatments
divdat <- data.frame(Cropping.system = substr(as.factor(rownames(divmat)), 1, 2),
                     Replicate = substr(as.factor(rownames(divmat)), 3, 3))
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(2,3,1,4)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Cropping.system.rep <- paste0(divdat$Cropping.system, divdat$Replicate)


# metadata
meta <- read.csv("Raw-data/study-sequences/Song-2018/SraRunTable_2.csv")
meta$Cropping.system.rep <- sapply(strsplit(meta$Library.Name, "-"), `[`, 1)
meta <- meta[,c("Cropping.system.rep", "SOM")]

divdat <- merge(divdat, meta, by="Cropping.system.rep")
divdat$ID <- paste0(divdat$Cropping.system.rep,".Fungus")

divdat$ID == rownames(divmat)





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

ENV$SOM2 <- ENV$SOM + rnorm(12,0.0001,0.0001)
vdat <- vegdist(t(OTU), method="jaccard")
mod<-dbrda(vdat ~ SOM +Condition(Replicate), data=ENV, na.action="na.omit")

h<-how(plots=Plots(strata=ENV$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=9999))
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/song_anova_fun.csv")
summod <- summary(mod)
write.csv(round(as.data.frame(summod[c("tot.chi", "partial.chi", "constr.chi", "unconst.chi")]),3), "Model-output/db-rda/song_chi_fun.csv")

bis <- as.data.frame(summod$biplot)
rownames(bis)
rownames(bis) <- c("SOC **")
bis$color <- c("black")


#  plot

cropgroups<-ENV$Cropping.system
levels(ENV$Cropping.system) <- c("Fallow soybean", "Continuous soybean", "Corn-soybean", "Wheat-soybean")
cols <- RColorBrewer::brewer.pal(n = 4, name = "Paired")

pdf("Figures/communities/Song_rda_fun.pdf", height=4, width=4)
par(mar=c(1.2,1.95,2,0.1))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=cols[cropgroups],
     pch=21, cex=1.25, xlim=c(-0.9, 1.5), 
     las=1, mgp = c(2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cropping system", 0.5,1.25, legend=levels(ENV$Cropping.system), col="black",
       pt.bg=cols, pch=21, pt.cex=1.25, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=bis$MDS1, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.25, 0, rownames(bis), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("B. Song ",italic("et al.")," (2018)")), side=3, adj=0, line=0.25, cex=1, font=4)

dev.off()






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
write.csv(mod.effects, "Model-output/Song-2018_permanova_fun.csv")
# examine treatment effects








