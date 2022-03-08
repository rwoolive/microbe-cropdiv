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



path1 <- "Raw-data/study-sequences/Cloutier-2020/seqs-deinterleaved" # path for fungi




################  FUNGI ################  
# asv taxonomy
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment_plant-pathogen.csv"))
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
plots <- 22:93
divmat <- t(taxabund[,plots])
colnames(divmat) <- taxabund$X.x
rownames(divmat) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat)
range(colSums(divmat)) # check that all asvs have at least 10 reads


# dataframe for diversity estimates and treatments
meta <- read.csv("Raw-data/study-sequences/Cloutier-2020/SraRunTable.csv")
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

write.csv(divdat, "Processed-data/Cloutier-diversity_pathogen.csv")


















