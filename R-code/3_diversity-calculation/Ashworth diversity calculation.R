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



path1 <- "Raw-data/study-sequences/Ashworth/All" # path for bacteria

authorM <- "Ashworth"
sampdateM <- c("year-1", "year-2") # sampling dates
locationM <- c("mtrec","recm") # sampling locations
crop.seqM <- c("cccc","cscs","ssss","ctctctct") # crop sequences
regionM <- "16S"




################  Bacteria ################  
# asv taxonomy
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqtab/3seqtab.csv"))

# species abundance curve
abunds <- abundances[,c(3:188)]
rownames(abunds) <- abundances$X.1
abunds2 <- t(abunds)

# mtrec curve
pdf("Figures/coverage/Ashworth-bac_mtrec.pdf", height=5, width=6)
rarecurve(abunds2[1:90,], step = 4000,  col = "blue", label = FALSE,  
          xlab = "Number of reads", ylab = "Number of ASVs")
dev.off()

# recm curve
pdf("Figures/coverage/Ashworth-bac_recm.pdf", height=5, width=6)
rarecurve(abunds2[91:186,], step = 4000,  col = "blue", label = FALSE,  
          xlab = "Number of reads", ylab = "Number of ASVs")
dev.off()

# Good's coverage
range(colSums(abunds2))
goodmat <- goods(abunds2)
mean(goodmat$goods)




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
                     ID = meta2$Sequence)
divdat$Cropping.system <- revalue(divdat$Cropping.system0, c("C>C>C>C"="Corn", "C>S>C>S"="Corn-soybean", "Ct>Ct>Ct>Ct"="Cotton", "S>S>S>S"="Soybean"))
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])


table(divdat$Cropping.system, divdat$Cover, divdat$Location, divdat$Yr)


# divmat_pa: matrix with ASVs as columns and plots as rows (proportional abundance)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_pa <- proportions(as.matrix(divmat),1)
range(rowSums(divmat_pa)) # check to see that all rows sum up to 1
divmat_pa <- data.frame(divmat_pa)


# divmat_presabs: matrix with ASVs as columns and plots as rows (presence/absence)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_presabs <- divmat
divmat_presabs[divmat_presabs>0] <- 1
range(rowSums(divmat_presabs)) # check to see that all rows sum up to integers
divmat[1:10,1:10] 
divmat_presabs[1:10,1:10] # compare to original divmat matrix



# look at summary of reads for each sample
summary(rowSums(divmat))
summary(rowSums(divmat_pa))
summary(rowSums(divmat_presabs))


########## Inverse Simpson's:
### (1/sum p_i^2)
invsimpson.div <- diversity(divmat_pa, "invsimpson") 
divdat$invsimpson.div <- invsimpson.div
hist(divdat$invsimpson.div)






######### RECM
site <- "_recm"
divdatr <- divdat[which(divdat$Location=="RECM"),]
divdatr$Replicate <- factor(divdatr$Replicate, levels(divdatr$Replicate)[c(1,2,3)]) # no fourth block at recm

write.csv(divdatr, "Processed-data/Ashworth_recm-diversity.csv")














######### MTREC
site <- "_MTREC"
divdatw <- divdat[which(divdat$Location=="MTREC"),]
divdatw$Cropping.system <- factor(divdatw$Cropping.system, levels(divdatw$Cropping.system)[c(1,2,4)]) # mtrec does not have cotton

write.csv(divdatw, "Processed-data/Ashworth_MTREC-diversity.csv")










######### Both RECM and MTREC, but exclude cotton
### ANOVA & plotting
site <- "_both-sites"
divdatb <- divdat[which(divdat$Cropping.system %in% c("Corn", "Soybean", "Corn-soybean")),]
divdatb$Cropping.system <- factor(divdatb$Cropping.system, levels(divdatb$Cropping.system)[c(1,2,4)]) # no cotton at mtrec

write.csv(divdatb, "Processed-data/Ashworth_both-sites-diversity.csv")





















