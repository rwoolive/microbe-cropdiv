
library(dada2); packageVersion("dada2")
library(ape)





##############################################
# Directory containing the fastq files after unzipping.
path <- "Raw-data/study-sequences/Gao/seqs-deinterleaved"



# Use the RDP Classifier
# Go to http://rdp.cme.msu.edu/classifier/classifier.jsp
# Citation:  Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.

# Under "Choose a gene:" select "16S rRNA training set 18"
# Then click the button that says "Choose file". Navigate to the 
# fasta file for ASVs with 10 reads or more. Then click "Submit".

# Save the output file in the taxonomy folder.
# Select a confidence threshold of 95%
# Open the file in excel, add headers, and save as csv

# Read in the file and trim the dataframe to the asvs that are 
# the correct kingdom classification 





taxa3 <- read.csv(paste0(path,"/taxonomy/fixrank_ASVs2_10readsormore.fa_classified.csv"))




# remove non-bacteria/archaea
table(taxa3$Kingdom) 
#taxa3[which(taxa3$Kingdom=="Eukaryota"),] # Genus: Zea
#nonbac <- which(taxa3$Kingdom=="Eukaryota")
#taxa4 <- taxa3[-nonbac,] # remove 
taxa4 <- taxa3

# remove the same asvs from the fasta file
bacfasta3 <- read.FASTA(paste0(path, "/seqtabs/ASVs2_10readsormore.fa"))
#bacfasta4 <- bacfasta3[-nonbac]
bacfasta4 <- bacfasta3
write.FASTA(bacfasta4, paste0(path, "/seqtabs/ASVs3_eukrem.fa"))

# remove the same asvs from the abundance file
abund <- read.csv(paste0(path, "/seqtabs/seqtab2_10readsormore.csv"))
#abund2 <- abund[-nonbac,]
abund2 <- abund
write.csv(abund2, paste0(path, "/seqtabs/seqtab3_eukrem.csv"))



length(which(is.na(taxa4$Kingdom)==TRUE)) # there are no unidentifiable asvs



write.csv(taxa4,  paste0(path, "/taxonomy/taxonomy-assignment.csv"))
taxa5 <- read.csv(paste0(path, "/taxonomy/taxonomy-assignment.csv"))

# number of asvs
write.csv(dim(taxa5)[1],paste0(path, "/taxonomy/asv-table_number-asvs.csv"))


# look at phyla represented
phyla <- as.data.frame(table(taxa5$Phylum))
phyla$percentage <- round((phyla$Freq/dim(taxa5)[1])*100,2)
write.csv(phyla, paste0(path, "/taxonomy/asv-table_phyla.csv"))

# look at classes represented
classes <- as.data.frame(table(taxa5$Class))
classes$percentage <- round((classes$Freq/dim(taxa5)[1])*100,2)
write.csv(classes, paste0(path, "/taxonomy/asv-table_classes.csv"))

# look at orders represented
orders <- as.data.frame(table(taxa5$Order))
orders$percentage <- round((orders$Freq/dim(taxa5)[1])*100,2)
write.csv(orders, paste0(path, "/taxonomy/asv-table_orders.csv"))

# look at families represented
family <- as.data.frame(table(taxa5$Family))
family$percentage <- round((family$Freq/dim(taxa5)[1])*100,2)
write.csv(family, paste0(path, "/taxonomy/asv-table_families.csv"))

# look at genera represented
genus <- as.data.frame(table(taxa5$Genus))
genus$percentage <- round((genus$Freq/dim(taxa5)[1])*100,2)
write.csv(genus, paste0(path, "/taxonomy/asv-table_genera.csv"))

# note that some of the species classifications look like genera, families
species <- table(taxa5$Species)
write.csv(species, paste0(path, "/taxonomy/asv-table_species.csv"))




########### subset to bacterial functional groups

taxa <- read.csv(paste0(path,"/taxonomy/taxonomy-assignment.csv"))
d <- data.frame(group=c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Acidobacteria", "Actinobacteria", "Firmicutes", "Planctomycetes", "Gemmatimonadetes", "Bacteroidetes", "Chloroflexi", "Verrucomicrobia", "Proteobacteria"),
                level=c("Class", "Class", "Class", "Phylum", "Phylum", "Phylum", "Phylum", "Phylum", "Phylum", "Phylum", "Phylum", "Phylum"))







# average abundance of bacterial groups across samples
abundances <- read.csv(paste0(path, "/seqtabs/seqtab3_eukrem.csv"))
sumabund <- colSums(abundances[,3:14])

# columns for plots in taxabund 
taxabund <- merge(taxa, abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE) 
plots <- 17:28

# create column in d for average relative abundance
d$avabund.percent <- rep(NA, dim(d)[1])

# merge taxa and abundance data into one df for groups 1-3 (orders)
for(i in 1:3){
  taxabund <- merge(taxa[taxa$Class==d$group[i],], abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE) 
  d$avabund.percent[i] <- round(mean(colSums(taxabund[,plots])/sumabund),3)*100
}

# merge taxa and abundance data into one df for groups 4-11 (phyla)
for(i in 4:12){
  taxabund <- merge(taxa[taxa$Phylum==d$group[i],], abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE) 
  d$avabund.percent[i] <-  round(mean(colSums(taxabund[,plots])/sumabund),3)*100
}

d$groupperc <- paste0(d$group, " (", d$avabund.percent, "%)")
d <- d[order(d$avabund.percent, decreasing = T),]
write.csv(d, paste0(path, "/taxonomy/bac-relative-abundances.csv"))














