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



path2 <- "Raw-data/study-sequences/Song-2018/seqs/fun" # path for fungi









################  Fungi ################  
# asv taxonomy
taxa <- read.csv(paste0(path2, "/taxonomy/taxonomy-assignment_amf-GLOM.csv"))
# asv abundances
abundances <- read.csv(paste0(path2, "/seqtabs/seqtab3_eukrem.csv"))


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abundances, 
                  by.x = "X.5", by.y = "X.1", 
                  all.x = TRUE, all.y = FALSE) # all.y = FALSE removes sequences in abundances that do not match up to the taxa dataset

# Kingdoms
table(taxabund$Kingdom)

# check ASV numbers by comparing X.3 and X.1 columns in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- 27:38
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

write.csv(divdat, "Processed-data/Song-diversity_fun_AMF.csv")




















