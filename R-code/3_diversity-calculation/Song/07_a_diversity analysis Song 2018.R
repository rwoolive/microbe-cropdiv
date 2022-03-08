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


# species abundance curve
abunds <- abundances[,c(3:14)]
rownames(abunds) <- abundances$X.1
abunds2 <- t(abunds)

pdf("Figures/coverage/Song-bac.pdf", height=5, width=6)
rarecurve(abunds2, step = 200,  col = "blue", label = FALSE,   
          xlab = "Number of reads", ylab = "Number of ASVs")
dev.off()

# Good's coverage
range(colSums(abunds2))
goodmat <- goods(abunds2)
mean(goodmat$goods)






# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by = "X.1", 
                  all.x = FALSE, all.y = TRUE) 

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
divdat$Cropping.system0 <- rep("Fallow soybean", dim(divdat)[1])
divdat$Cropping.system0[which(divdat$Cropping.system=="SC")] <- "Continuous soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="WS")] <- "Wheat-soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="CS")] <- "Corn-soybean"
divdat$Cropping.system0 <- as.factor(divdat$Cropping.system0)
divdat$Cropping.system0 <- factor(divdat$Cropping.system0, levels(divdat$Cropping.system0)[c(3,1,2,4)])
divdat$Cropping.system <- divdat$Cropping.system0
divdat$Replicate <- as.factor(divdat$Replicate)

table(divdat$Cropping.system)





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

write.csv(divdat, "Processed-data/Song-diversity_bac.csv")



  












################  Fungi ################  
# asv taxonomy
taxa <- read.csv(paste0(path2, "/taxonomy/taxonomy-assignment.csv"))
# asv abundances
abundances <- read.csv(paste0(path2, "/seqtabs/seqtab3_eukrem.csv"))



# species abundance curve
abunds <- abundances[,c(3:14)]
rownames(abunds) <- abundances$X.1
abunds2 <- t(abunds)

pdf("Figures/coverage/Song-fun.pdf", height=5, width=6)
rarecurve(abunds2, step = 200,  col = "blue", label = FALSE,  
          xlab = "Number of reads", ylab = "Number of ASVs")
dev.off()

# Good's coverage
range(colSums(abunds2))
goodmat <- goods(abunds2)
mean(goodmat$goods)








# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.5", by.y = "X.1", 
                  all.x = TRUE, all.y = FALSE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.3 and X.1 columns in taxabund


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
divdat$Cropping.system0 <- rep("Fallow soybean", dim(divdat)[1])
divdat$Cropping.system0[which(divdat$Cropping.system=="SC")] <- "Continuous soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="WS")] <- "Wheat-soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="CS")] <- "Corn-soybean"
divdat$Cropping.system0 <- as.factor(divdat$Cropping.system0)
divdat$Cropping.system0 <- factor(divdat$Cropping.system0, levels(divdat$Cropping.system0)[c(3,1,2,4)])
divdat$Cropping.system <- divdat$Cropping.system0
divdat$Replicate <- as.factor(divdat$Replicate)





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

write.csv(divdat, "Processed-data/Song-diversity_fun.csv")




















