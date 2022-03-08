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



path1 <- "Raw-data/study-sequences/Gao/seqs-deinterleaved" # path for bacteria




################  Bacteria ################  
# asv taxonomy
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqtabs/seqtab3_eukrem.csv"))


# species abundance curve
abunds <- abundances[,c(3:14)]
rownames(abunds) <- abundances$X.1
colnames(abunds) <- divdat$ID2 # run meta code below to obtain this
abunds2 <- t(abunds)

pdf("Figures/coverage/Gao-bac.pdf", height=5, width=6)
rarecurve(abunds2, step = 20,  col = "blue", label = FALSE,  
          xlab = "Number of reads", ylab = "Number of ASVs")
dev.off()
# Good's coverage
range(colSums(abunds2))
goodmat <- goods(abunds2)
mean(goodmat$goods)




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
meta <- read.csv("Raw-data/study-sequences/Gao/SraRunTable_updated.csv")
meta2 <- meta[which(meta$Experiment %in% rownames(divmat)),]
rownames(meta2) <- meta2$Experiment
meta2 <- meta2[rownames(divmat),]
meta2$Sample.Name2 <- sapply(strsplit(meta2$Sample.Name, " "), `[`, 1)
meta2$Sample.Name3 <- rep("Wheat", dim(meta2)[1])
meta2$Sample.Name3[which(meta2$Sample.Name2=="5-yearAgroforestry-cropalley")] <- "AF5"
meta2$Sample.Name3[which(meta2$Sample.Name2=="9-yearAgroforestry-cropalley")] <- "AF9"
meta2$Sample.Name3[which(meta2$Sample.Name2=="14-yearAgroforestry-cropalley")] <- "AF14"
meta2$Sample.Name3 <- as.factor(meta2$Sample.Name3)
meta2$Sample.Name3 <- factor(meta2$Sample.Name3, levels(meta2$Sample.Name3)[c(4,2,3,1)])
meta2$Rep <- c(2,1,1,3,3,2,1,3,2,1,3,2)
divdat <- data.frame(Cropping.system = meta2$Sample.Name3,
                     Replicate = meta2$Rep, 
                     ID = meta2$Experiment)
divdat$Cropping.system
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$ID2 <- paste(divdat$Cropping.system, divdat$Replicate, sep="_")

table(divdat$Cropping.system)

# check that datasets are ordered correctly
divdat$ID
row.names(divmat)



# divmat_pa: matrix with ASVs as columns and plots as rows (proportional abundance)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_pa <- proportions(as.matrix(divmat),1)
range(rowSums(divmat_pa)) # check to see that all rows sum up to 1
divmat_pa <- data.frame(divmat_pa)




# look at summary of reads for each sample
summary(rowSums(divmat))
summary(rowSums(divmat_pa))




########## Inverse Simpson's:
### (1/sum p_i^2)
invsimpson.div <- diversity(divmat_pa, "invsimpson") 
divdat$invsimpson.div <- invsimpson.div
hist(divdat$invsimpson.div)




write.csv(divdat, "Processed-data/Gao-diversity.csv")


















