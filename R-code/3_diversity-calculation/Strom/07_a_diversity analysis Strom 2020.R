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
taxa <- read.csv(paste0(path1, "/taxonomy/taxonomy-assignment-with-ft.csv"))
# asv abundances
abundances <- read.csv(paste0(path1, "/seqtabs/seqtab3_eukrem.csv"))


# species abundance curve
abunds <- abundances[,c(3:98)]
rownames(abunds) <- abundances$X.1
abunds <- abunds[,order(as.numeric(substr(colnames(abunds), start=4,stop=10)))]
divdat2 <- divdat[order(as.numeric(substr(divdat$ID, start=4,stop=10))),]
colnames(abunds) <- divdat2$ID2 # run meta code below to obtain this
abunds2 <- t(abunds)


pdf("Figures/coverage/Strom-fun.pdf", height=5, width=6)
rarecurve(abunds2, step = 200,  col = "blue", label = FALSE,  
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
                     Year = meta2$year,
                     Season = meta2$season,
                     Replicate = meta2$block, 
                     ID = meta2$Experiment)
divdat$Cropping.system0 <- rep("Corn", dim(divdat)[1])
divdat$Cropping.system0[which(divdat$Cropping.system=="Ca")] <- "(Corn)-soybean"
divdat$Cropping.system0[which(divdat$Cropping.system=="Sa")] <- "Corn-(soybean)"
divdat$Cropping.system0[which(divdat$Cropping.system=="Ss")] <- "Soybean"
divdat$Cropping.system0 <- as.factor(divdat$Cropping.system0)
divdat$Cropping.system0 <- factor(divdat$Cropping.system0, levels(divdat$Cropping.system0)[c(2,1,4,3)])
divdat$Cropping.system <- divdat$Cropping.system0
divdat$Year <- as.factor(divdat$Year)
divdat$Season <- as.factor(divdat$Season)
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(3,2,1)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$ID2 <- paste0(divdat$Cropping.system, "_rep-", divdat$Replicate, "_", divdat$Year, "_", divdat$Season)

table(divdat$Cropping.system, divdat$Year, divdat$Season)




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



write.csv(divdat, "Processed-data/Strom-diversity.csv")


















