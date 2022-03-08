


##############################################
# DADA2: https://benjjneb.github.io/dada2/index.html



##############################################
# Starting point
# Raw sequence data in the form of compressed .fastq files 
# These files were obtained directly from the authors (J. DeBruyn)
# Place files into Raw-data/study-sequences directory



##############################################
# install DADA2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")
# Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)



##############################################
# load libraries
library(dada2)
packageVersion("dada2") # version 1.16
library(seqinr)
library(ggplot2)




##############################################
# Sequences from separate locations, years, and cropping sequences were filtered  
# and trimmed in previous R codes. 
# In the following code we will merge sequences into one file then complete
# sequence processing according to the DADA tutorial.







###################
# Bacteria
###################

pathM <- "Raw-data/study-sequences/"
authorM <- "Ashworth"
sampdateM <- c("year-1", "year-2") # sampling dates
locationM <- c("mtrec","recm") # sampling locations
crop.seqM <- c("cccc","cscs","ssss","ctctctct") # crop sequences
regionM <- "16S"


# read in metadata 
metdat <- read.csv(paste0(pathM,authorM,"/Ash-Metadata.csv"))

### read in sequence tables from each sampling period and clean up sample names
# mtrec
# old way of loading data: load(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[1],"-",crop.seqM[2],".RData"))
# use cgwtools::lsdata() to get object names
seqtab.nochim_M1cccc <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[1],"-",crop.seqM[1],".RData"))[["seqtab.nochim_M1cccc"]]
seqtab.nochim_M1cccc <- t(seqtab.nochim_M1cccc)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M1cccc))
rownames(seqtab.nochim_M1cccc) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_M1cscs <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[1],"-",crop.seqM[2],".RData"))[["seqtab.nochim_M2cscs"]]
seqtab.nochim_M1cscs <- t(seqtab.nochim_M1cscs)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M1cscs))
rownames(seqtab.nochim_M1cscs) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_M1ssss <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[1],"-",crop.seqM[3],".RData"))[["seqtab.nochim_M1ssss"]]
seqtab.nochim_M1ssss <- t(seqtab.nochim_M1ssss)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M1ssss))
rownames(seqtab.nochim_M1ssss) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_M2cccc <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[2],"-",crop.seqM[1],".RData"))[["seqtab.nochim_M2cccc"]]
seqtab.nochim_M2cccc <- t(seqtab.nochim_M2cccc)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M2cccc))
rownames(seqtab.nochim_M2cccc) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_M2cscs <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[2],"-",crop.seqM[2],".RData"))[["seqtab.nochim_M2cscs"]]
seqtab.nochim_M2cscs <- t(seqtab.nochim_M2cscs)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M2cscs))
rownames(seqtab.nochim_M2cscs) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_M2ssss <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[1],"-",sampdateM[2],"-",crop.seqM[3],".RData"))[["seqtab.nochim_M2ssss"]]
seqtab.nochim_M2ssss <- t(seqtab.nochim_M2ssss)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_M2ssss))
rownames(seqtab.nochim_M2ssss) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])


# recm
seqtab.nochim_R1cccc <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[1],"-",crop.seqM[1],".RData"))[["seqtab.nochim_R1cccc"]]
seqtab.nochim_R1cccc <- t(seqtab.nochim_R1cccc)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R1cccc))
rownames(seqtab.nochim_R1cccc) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R1cscs <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[1],"-",crop.seqM[2],".RData"))[["seqtab.nochim_R1cscs"]]
seqtab.nochim_R1cscs <- t(seqtab.nochim_R1cscs)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R1cscs))
rownames(seqtab.nochim_R1cscs) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R1ssss <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[1],"-",crop.seqM[3],".RData"))[["seqtab.nochim_R1ssss"]]
seqtab.nochim_R1ssss <- t(seqtab.nochim_R1ssss)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R1ssss))
rownames(seqtab.nochim_R1ssss) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R1ctctctct <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[1],"-",crop.seqM[4],".RData"))[["seqtab.nochim_R2ctctctct"]]
seqtab.nochim_R1ctctctct <- t(seqtab.nochim_R1ctctctct)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R1ctctctct))
rownames(seqtab.nochim_R1ctctctct) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R2cccc <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[2],"-",crop.seqM[1],".RData"))[["seqtab.nochim_R2cccc"]]
seqtab.nochim_R2cccc <- t(seqtab.nochim_R2cccc)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R2cccc))
rownames(seqtab.nochim_R2cccc) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R2cscs <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[2],"-",crop.seqM[2],".RData"))[["seqtab.nochim_R2cscs"]]
seqtab.nochim_R2cscs <- t(seqtab.nochim_R2cscs)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R2cscs))
rownames(seqtab.nochim_R2cscs) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R2ssss <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[2],"-",crop.seqM[3],".RData"))[["seqtab.nochim_R2ssss"]]
seqtab.nochim_R2ssss <- t(seqtab.nochim_R2ssss)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R2ssss))
rownames(seqtab.nochim_R2ssss) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])

seqtab.nochim_R2ctctctct <- loadToEnv(paste0(pathM,authorM,"/ash fastq fileslist-",locationM[2],"-",sampdateM[2],"-",crop.seqM[4],".RData"))[["seqtab.nochim_R2ctctctct"]]
seqtab.nochim_R2ctctctct <- t(seqtab.nochim_R2ctctctct)
ourrows <- which(metdat$Sequence %in% rownames(seqtab.nochim_R2ctctctct))
rownames(seqtab.nochim_R2ctctctct) <- paste0(regionM, "-", metdat$Location[ourrows], "-", metdat$SequenceLet[ourrows], "-", metdat$Cover[ourrows], "-Rep", metdat$Rep[ourrows], "-Yr", metdat$Yr[ourrows])



# The sequence table is a matrix with rows corresponding to (and named by) the 
# samples, and columns corresponding to (and named by) the sequence variants. 




# inspect sequence lengths:
table(nchar(colnames(seqtab.nochim_M1cccc)))
table(nchar(colnames(seqtab.nochim_M1cscs)))
table(nchar(colnames(seqtab.nochim_M1ssss)))
table(nchar(colnames(seqtab.nochim_M2cscs)))
table(nchar(colnames(seqtab.nochim_M2ssss)))
table(nchar(colnames(seqtab.nochim_R1cccc)))
table(nchar(colnames(seqtab.nochim_R1cscs)))
table(nchar(colnames(seqtab.nochim_R1ssss)))
table(nchar(colnames(seqtab.nochim_R1ctctctct)))
table(nchar(colnames(seqtab.nochim_R2cccc)))
table(nchar(colnames(seqtab.nochim_R2cscs)))
table(nchar(colnames(seqtab.nochim_R2ssss)))
table(nchar(colnames(seqtab.nochim_R2ctctctct)))

# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.

# NOTE: Sequences that are much longer or shorter than expected may be the 
# result of non-specific priming. You can remove non-target-length sequences 
# from your sequence table: 
# eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 
# This is analogous to “cutting a band” in-silico to get amplicons of the 
# targeted length.


# Combine datasets from different sampling periods:
seqtab.nochim_combo <- mergeSequenceTables(seqtab.nochim_M1cccc, seqtab.nochim_M1cscs, seqtab.nochim_M1ssss, 
                                    seqtab.nochim_M2cccc, seqtab.nochim_M2cscs, seqtab.nochim_M2ssss,
                                    seqtab.nochim_R1cccc, seqtab.nochim_R1cscs, seqtab.nochim_R1ssss, seqtab.nochim_R1ctctctct,
                                    seqtab.nochim_R2cccc, seqtab.nochim_R2cscs, seqtab.nochim_R2ssss, seqtab.nochim_R2ctctctct)


##############################################
# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made 
# it through each step in the pipeline:

getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(#out, 
               #sapply(dadaFs, getN), 
               #sapply(dadaRs, getN), 
               #sapply(mergers, getN), 
               rowSums(seqtab.nochim_combo))#,
               #round(rowSums(seqtab.nochim)/out[,1]*100, 1))
colnames(track) <- c(#"input", 
                     #"filtered", 
                     #"denoisedF", 
                     #"denoisedR", 
                     #"merged", 
                     "nonchim")#,
                     #"final_perc_reads_retained")
head(track)
track <- as.data.frame(track)
write.csv(track, paste0("Raw-data/study-sequences/",authorM,"/All/track.csv"))


track <- read.csv(paste0("Raw-data/study-sequences/",authorM,"/All/track.csv"))
sum(track$nonchim[1:90]) # number of reads from MTREC: 7387600
sum(track$nonchim[91:186]) # number of reads from RECM: 7513455
sum(track$nonchim) # total number of reads: 14901055


# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)


# Looks good! We kept the majority of our raw reads, and there is no over-large 
# drop associated with any single step.

# NOTE: This is a great place to do a last sanity check. Outside of filtering, 
# there should no step in which a majority of reads are lost. If a majority of 
# reads failed to merge, you may need to revisit the truncLen parameter used in 
# the filtering step and make sure that the truncated reads span your amplicon. 
# If a majority of reads were removed as chimeric, you may need to revisit the 
# removal of primers, as the ambiguous nucleotides in unremoved primers interfere 
# with chimera identification.


##############################################
# Save the seqtab file with no chimeras
seqtab.nochim_combo <- t(seqtab.nochim_combo)
write.csv(seqtab.nochim_combo, file=paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/0seqtab.csv"))




##############################################
# Remove sequences that are found within the PCR blank,
# and then remove sequences with less than 10 reads

seqtab.nochim_combo <- t(seqtab.nochim_combo)
seqtab.nochim_combo_df <- as.data.frame(seqtab.nochim_combo)
# determine if any asvs were detected in pcr blank
#unique(seqtab.nochim_combo_df$`PCRblank-16S_2020-09`); which(seqtab.nochim_df$`PCRblank-16S_2020-09`>0) 
# there was one asv that was detected in the blank, but it had 2 reads so it  
# will be removed anyway

# remove column with pcr blank
#blankcol <- which(colnames(seqtab.nochim_df)=="PCRblank-16S_2020-09")
seqtab.nochim_df2 <- seqtab.nochim_combo_df#[,-blankcol]


#### remove asvs with less than 10 reads across samples
seqtab.nochim_df2$sum <- rowSums(seqtab.nochim_df2)
length(seqtab.nochim_df2$sum) # check against 1seqtab-dim.csv

# number of asvs with less than 10 reads
write.csv(length(which(seqtab.nochim_df2$sum<10)), file=paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/1lessthan10reads.csv"))
# number of asvs with > 10 reads
write.csv(length(which(seqtab.nochim_df2$sum>=10)), file=paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/1atleast10reads.csv"))
# proportion of asvs retained: 0.3468469
write.csv(length(which(seqtab.nochim_df2$sum>=10))/length(seqtab.nochim_df2$sum), file=paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/1proportionretained.csv"))
 
# remove
seqtab.nochim_df3 <- seqtab.nochim_df2[-which(seqtab.nochim_df2$sum<10),]
hist(seqtab.nochim_df3$sum); range(seqtab.nochim_df3$sum) # check

write.csv(seqtab.nochim_df3, file=paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/2seqtab.csv"))



##############################################
# Summary tables
seqtab.nochim4 <- t(seqtab.nochim_df3)
asv_seqs <- colnames(seqtab.nochim4)
asv_headers <- vector(dim(seqtab.nochim4)[2], mode="character")

for (i in 1:dim(seqtab.nochim4)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/2seqtab.fa"))

# count table:
asv_tab <- t(seqtab.nochim4)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0("Raw-data/study-sequences/",authorM,"/All/seqtab/2counts.tsv"), sep="\t", quote=F, col.names=NA)






