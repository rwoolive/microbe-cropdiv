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
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment_all-saprotroph.csv"))
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
plots <- 22:117
divmat <- t(taxabund[,plots])
colnames(divmat) <- taxabund$X.x
rownames(divmat) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat)
range(colSums(divmat)) # check that all asvs have at least 10 reads


# dataframe for diversity estimates and treatments
meta <- read.csv("Raw-data/study-sequences/Strom-2020/SraRunTable.csv")
meta2 <- meta[which(meta$Experiment %in% rownames(divmat)),]
divdat <- data.frame(Cropping.system = meta2$crop_rotation,
                     Replicate = meta2$block, 
                     ID = meta2$Experiment,
                     Season = meta2$collection_date,
                     Year = meta2$collection_date)
divdat$Cropping.system0 <- rep("Corn", dim(divdat)[1])
divdat$Cropping.system0[which(divdat$Cropping.system=="Ca")] <- "(Corn)-soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="Sa")] <- "Corn-(soybean)"
divdat$Cropping.system0[which(divdat$Cropping.system=="Ss")] <- "Soybean"
divdat$Cropping.system0 <- as.factor(divdat$Cropping.system0)
divdat$Cropping.system0 <- factor(divdat$Cropping.system0, levels(divdat$Cropping.system0)[c(2,1,4,3)])
divdat$Cropping.system <- divdat$Cropping.system0
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$ID2 <- paste(divdat$Cropping.system, divdat$Replicate, sep="_")
divdat$Season0 <- rep("Spring", dim(divdat)[1])
divdat$Season0[which(substr(divdat$Season, 6,7) %in% c("07"))] <- "Summer"
divdat$Season0[which(substr(divdat$Season, 6,7) %in% c("10"))] <- "Fall"
divdat$Season0 <- as.factor(divdat$Season0)
divdat$Season0 <- factor(divdat$Season0, levels(divdat$Season0)[c(2,3,1)])
divdat$Season <- divdat$Season0
divdat$Year0 <- rep(2015, dim(divdat)[1])
divdat$Year0[which(substr(divdat$Year, 1,4) %in% c("2016"))] <- 2016
divdat$Year0 <- as.factor(divdat$Year0)
divdat$Year <- divdat$Year0

table(divdat$Cropping.system)
table(divdat$Season)
table(divdat$Year)




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



write.csv(divdat, "Processed-data/Strom-diversity_all-saprotroph.csv")


















