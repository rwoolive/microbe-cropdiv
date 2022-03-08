
library(dada2); packageVersion("dada2")
library(ape)





##############################################
# Directory containing the fastq files after unzipping.
path <- "Raw-data/study-sequences/Cloutier-2020/seqs-deinterleaved"



################# Taxonomic assignment: Bacteria ################# 

### read in output files from DADA2 sequence processing code
### these files have the asv sequences in the first column and the
### sequence abundances for each sample in the columns thereafter
# fungi
fun <- read.csv(paste0(path,"/seqtabs/seqtab1.csv"), row.names=1)


##### determine if any asvs were detected in pcr blank
# no pcr blank in this dataset
#unique(fun$PCRblank.16S); which(fun$PCRblank.16S>0) 


##### remove pcr blank
#fun <- fun[,-81]




# Use the RDP Classifier
# Go to http://rdp.cme.msu.edu/classifier/classifier.jsp
# Citation:  Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.

# Under "Choose a gene:" select "16S rRNA training set 18" for bacteria
# or "UNITE Fungal ITS trainset 07-04-2014" for fungi
# Then click the button that says "Choose file". Navigate to the 
# fasta file for ASVs with 10 reads or more. Then click "Submit".

# Save the output file in the taxonomy folder.
# Select a confidence threshold of 95%
# Open the file in excel, add headers, and save as csv

# Read in the file and trim the dataframe to the asvs that are 
# the correct kingdom classification 



funtax <- read.csv(paste0(path,"/taxonomy/fixrank_ASVs2_10readsormore.fa_classified.csv"))





# fungi
table(funtax$Kingdom)  # all fungi
#funtax[which(funtax$Kingdom=="Eukaryota"),] # Genus: Zea
#nonfun <- which(funtax$Kingdom=="Eukaryota")
funtax2 <- funtax#[-nonfun,] # remove 

# remove the same asvs from the fasta file
funfasta3 <- read.FASTA(paste0(path, "/seqtabs/ASVs2_10readsormore.fa"))
funfasta4 <- funfasta3#[-nonfun]
write.FASTA(funfasta4, paste0(path, "/seqtabs/ASVs3_eukrem.fa"))

# remove the same asvs from the abundance file
abund <- read.csv(paste0(path, "/seqtabs/seqtab2_10readsormore.csv"))
abund2 <- abund#[-nonfun,]
write.csv(abund2, paste0(path, "/seqtabs/seqtab3_eukrem.csv"))

length(which(is.na(funtax2$Kingdom)==TRUE)) # there are no unidentifiable asvs

write.csv(funtax2,  paste0(path, "/taxonomy/taxonomy-assignment.csv"))
funtax3 <- read.csv(paste0(path, "/taxonomy/taxonomy-assignment.csv"))

# number of asvs
write.csv(dim(funtax3)[1],paste0(path, "/taxonomy/asv-table_number-asvs.csv"))




### biological groups represented





# fungi

# look at phyla represented
phyla <- as.data.frame(table(funtax3$Phylum))
phyla$percentage <- round((phyla$Freq/dim(funtax3)[1])*100,2)
write.csv(phyla, paste0(path, "/taxonomy/asv-table_phyla.csv"))

# look at classes represented
classes <- as.data.frame(table(funtax3$Class))
classes$percentage <- round((classes$Freq/dim(funtax3)[1])*100,2)
write.csv(classes, paste0(path, "/taxonomy/asv-table_classes.csv"))

# look at orders represented
orders <- as.data.frame(table(funtax3$Order))
orders$percentage <- round((orders$Freq/dim(funtax3)[1])*100,2)
write.csv(orders, paste0(path, "/taxonomy/asv-table_orders.csv"))

# look at families represented
family <- as.data.frame(table(funtax3$Family))
family$percentage <- round((family$Freq/dim(funtax3)[1])*100,2)
write.csv(family, paste0(path, "/taxonomy/asv-table_families.csv"))

# look at genera represented
genus <- as.data.frame(table(funtax3$Genus))
genus$percentage <- round((genus$Freq/dim(funtax3)[1])*100,2)
write.csv(genus, paste0(path, "/taxonomy/asv-table_genera.csv"))

# note that some of the species classifications look like genera, families
species <- as.data.frame(table(funtax3$Species))
species$percentage <- round((species$Freq/dim(funtax3)[1])*100,2)
write.csv(species, paste0(path, "/taxonomy/asv-table_species.csv"))















########### Assign function using Functional trait database
# functional trait data from Polme et al. 2021 Fungal Diversity
# https://link.springer.com/article/10.1007/s13225-020-00466-2


ftrait <- read.csv("Raw-data/Fungal-trait-database/data.csv")

taxa <- read.csv(paste0(path, "/taxonomy/taxonomy-assignment.csv"))
taxa$Primary_lifestyle <- rep(NA, dim(taxa)[1])
taxa$Secondary_lifestyle <- rep(NA, dim(taxa)[1])

for(i in 1:dim(taxa)[1]){
  if(taxa$Genus[i] %in% ftrait$GENUS){
    taxa$Primary_lifestyle[i] <- ftrait$primary_lifestyle[which(ftrait$GENUS==taxa$Genus[i])]
    taxa$Secondary_lifestyle[i] <- ftrait$Secondary_lifestyle[which(ftrait$GENUS==taxa$Genus[i])]
  }
}


# functional trait distribution
ftraittab <- as.data.frame(table(taxa$Primary_lifestyle))
write.csv(ftraittab, paste0(path, "/taxonomy/5a_primary-lifestyle-asvs.csv"))

# create data subsets for certain functional groups
write.csv(taxa[which(taxa$Primary_lifestyle=="arbuscular_mycorrhizal"),], paste0(path, "/taxonomy/taxonomy-assignment_amf.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="soil_saprotroph"),], paste0(path, "/taxonomy/taxonomy-assignment_soil-saprotroph.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="litter_saprotroph"),], paste0(path, "/taxonomy/taxonomy-assignment_litter-saprotroph.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle %in% c("soil_saprotroph", "litter_saprotroph", "dung_saprotroph", "nectar/tap_saprotroph", "pollen_saprotroph", "unspecified_saprotroph", "wood_saprotroph")),], paste0(path, "/taxonomy/taxonomy-assignment_all-saprotroph.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="plant_pathogen"),], paste0(path, "/taxonomy/taxonomy-assignment_plant-pathogen.csv"))


write.csv(taxa, paste0(path,"/taxonomy/taxonomy-assignment-with-ft.csv"))




#### Subset to glomeromycota to explore AMF


taxa5 <- read.csv(paste0(path,"/taxonomy/taxonomy-assignment-with-ft.csv"))
as.data.frame(table(taxa5$Primary_lifestyle))
amftax <- taxa5[which(taxa5$Phylum=="Glomeromycota"),]

write.csv(amftax, paste0(path,"/taxonomy/taxonomy-assignment_amf-GLOM.csv"))








########### subset to fungal functional groups and determine average relative abundances

taxa <- read.csv(paste0(path,"/taxonomy/taxonomy-assignment-with-ft.csv"))
d <- data.frame(group=c("arbuscular_mycorrhizal", "soil_saprotroph", "litter_saprotroph", "dung_saprotroph", "nectar/tap_saprotroph", "pollen_saprotroph", "unspecified_saprotroph", "wood_saprotroph", "plant_pathogen", "saprotroph",
                        "Ascomycota", "Basidiomycota", "Zygomycota", "Chytridiomycota", "Glomeromycota"),
                level=c(rep("Primary_lifestyle",10), rep("Phylum",5)))

# average abundance of fungal groups across samples
abundances <- read.csv(paste0(path, "/seqtabs/seqtab3_eukrem.csv"))
sumabund <- colSums(abundances[,3:74])

# columns for plots in taxabund
taxabund <- merge(taxa, abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE) 
plots <- 22:93

# create column in d for average relative abundance
d$avabund.percent <- rep(NA, dim(d)[1])

# merge taxa and abundance data into one df for groups 1-10 (primary lifestyle)
for(i in 1:10){
  taxabund <- merge(taxa[which(taxa$Primary_lifestyle==d$group[i]),], abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE)
  d$avabund.percent[i] <- round(mean(colSums(taxabund[,plots])/sumabund),3)*100
}

# merge taxa and abundance data into one df for groups 4-11 (phyla)
for(i in 11:15){
  taxabund <- merge(taxa[taxa$Phylum==d$group[i],], abundances, by.x = "X.1", by.y = "X.1", all.x = TRUE, all.y = FALSE) 
  d$avabund.percent[i] <-  round(mean(colSums(taxabund[,plots])/sumabund),3)*100
}

d$avabund.percent[which(d$group=="saprotroph")] <- sum(d$avabund.percent[which(d$group %in% c("soil_saprotroph", "litter_saprotroph", "dung_saprotroph", "nectar/tap_saprotroph", "pollen_saprotroph", "unspecified_saprotroph", "wood_saprotroph"))])

d$groupperc <- paste0(d$group, " (", d$avabund.percent, "%)")
d <- d[order(d$avabund.percent, decreasing = T),]
write.csv(d, paste0(path, "/taxonomy/fun-relative-abundances.csv"))




