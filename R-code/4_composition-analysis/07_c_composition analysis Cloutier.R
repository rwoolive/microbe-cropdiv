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







################  Fungi ################  
path1 <- "Raw-data/study-sequences/Cloutier-2020" # path for bacteria

author <- "Cloutier"
group <- "fun"
region <- "ITS"




# asv taxonomy
taxa <- read.csv(paste0(path1, "/seqs-deinterleaved/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqs-deinterleaved/seqtabs/seqtab3_eukrem.csv"))


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.1", by.y = "X.1", 
                  all.x = TRUE, all.y = TRUE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.1 and X.x columns in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- 19:90 # columns that contain abundances
divmat <- taxabund[,plots]
divmat2 <- t(divmat)
colnames(divmat2) <- taxabund$X.x
rownames(divmat2) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat2)
range(colSums(divmat)) # check that all asvs have at least 10 reads


# dataframe for diversity estimates and treatments
meta <- read.csv(paste0(path1, "/SraRunTable.csv"))
meta2 <- meta[which(meta$Experiment %in% rownames(divmat)),]

divdat <- data.frame(Cover = as.factor(sapply(strsplit(basename(meta2$Library.Name), " CC Block"), `[`, 1)),
                     Replicate = sapply(strsplit(basename(meta2$Library.Name), " CC Block "), `[`, 2), 
                     Season = sapply(strsplit(basename(meta2$Library.Name), " CC Block "), `[`, 2),
                     ID = meta2$Experiment)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(5,3,4,6,7,8,9,1,2)])
divdat$Replicate <- as.factor(sapply(strsplit(basename(divdat$Rep), " "), `[`, 1))
divdat$Season <- as.factor(sapply(strsplit(basename(divdat$Season), " "), `[`, 2))
table(divdat$Cover, divdat$Season) # replication
divdat$ID2 <- paste(divdat$Cover, divdat$Replicate, divdat$Season, sep="_")
divdat <- divdat[order(divdat$ID),]

table(divdat$Cover, divdat$Season)


# merge with carbon data
cdat <- read.csv(paste0(path1, "/CoverCropFungal-CarbonData.csv"))
divdat$ID2[1]
cdat$X[1]
cdat$Cover <- sapply(strsplit(cdat$X, "-"), `[`, 1)
cdat$Cover[which(cdat$Cover=="3SppN")] <- "3 Spp mix"
cdat$Cover[which(cdat$Cover=="6Spp")] <- "6 Spp mix"
cdat$Replicate <- sapply(strsplit(cdat$X, "-"), `[`, 2)
cdat$Season <- sapply(strsplit(cdat$X, "-"), `[`, 3)
cdat$Season[which(cdat$Season=="1")] <- "Spring"
cdat$Season[which(cdat$Season=="2")] <- "Summer"
cdat$ID2 <- paste(cdat$Cover, cdat$Replicate, cdat$Season, sep="_")
cdat <- cdat[,c("ID2", "OM", "POXC")]
divdat$ID2[1]
cdat$ID2[1]

divdat <- merge(divdat, cdat, by = "ID2") 

divdat <- divdat[order(divdat$ID),]

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
mod<-dbrda(vdat ~ OM + POXC + Condition(Replicate), data=ENV, na.action="na.omit")

h<-how(plots=Plots(strata=ENV$Replicate))

modout <- anova(mod, by="margin", permuations=how(nperm=999))
write.csv(round(as.data.frame(modout),3), "Model-output/db-rda/cloutier_anova.csv")
summod <- summary(mod)
write.csv(round(as.data.frame(summod[c("tot.chi", "partial.chi", "constr.chi", "unconst.chi")]),3), "Model-output/db-rda/cloutier_chi.csv")

bis <- as.data.frame(summod$biplot)
rownames(bis)
rownames(bis) <- c("SOC", "POXC")
bis$color <- c("black","black")


#  plot

covgroups<-ENV$Cover
ssngroups <-ENV$Season
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set3")


pdf("Figures/communities/Cloutier_rda.pdf", height=4, width=4)
par(mar=c(1.2,1.2,2,0.2))
plot(scores(mod)$sites[,1],scores(mod)$sites[,2],
     bg=cols[covgroups],
     pch=c(21,22)[ssngroups], cex=1.25, xlim=c(-2, 3.5), 
     las=1, mgp = c(2.2, 0.3, 0), tcl=0.2)
abline(v=0, lty=3, lwd=0.5); abline(h=0, lty=3, lwd=0.5)
legend(title="Cover", 2.2,-0.1, legend=levels(ENV$Cover), col="black",
       pt.bg=cols, pch=21, pt.cex=1.25, cex=0.75)
legend(title="Season", 2,1.8, legend=levels(ENV$Season), 
       pch=c(21,22), pt.cex=1, cex=0.75)
arrows(x0=0, x1=bis$dbRDA1, y0=0, y1=bis$dbRDA2, col=bis$color, lwd=2, length=0.1)
text(bis$dbRDA1+sign(bis$dbRDA1)*0.3, bis$dbRDA2+sign(bis$dbRDA2)*0.15, rownames(bis), col=bis$color, font=2, cex=0.75)
mtext(expression(paste("A. Cloutier ",italic("et al.")," (2020)")), side=3, adj=0, line=0.25, cex=1, font=4)
dev.off()
















#######################################


### permanova 
colnames(divmat_pa)[1:5]; rownames(divmat_pa)[1:5] # asvs as columns, samples as rows
colnames(divdat)[1:5]; rownames(divdat)[1:5] # treatment levels as columns, samples as rows


mod.effects <- data.frame(effect=c("Cover", "Season", 
                                   "Cover:Season",  
                                   "Residual", "Total")) #  effects


# permanova model
permmod <- adonis(divmat_pa ~ Cover*Season, 
                  data=divdat, 
                  method="bray", 
                  strata=divdat$Replicate) 
permmod2 <- round(as.data.frame(permmod$aov.tab), 2)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$result <- ssp



write.csv(mod.effects, paste("Model-output/",author,"_permanova_",group, region,".csv"))
# examine for treatment effects




