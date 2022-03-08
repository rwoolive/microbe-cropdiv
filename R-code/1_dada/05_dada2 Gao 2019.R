


##############################################
# DADA2: https://benjjneb.github.io/dada2/index.html



##############################################
# Starting point
# Raw sequence data in the form of compressed .fastq files 
# These files were obtained from National Center for Biotechnology Information under the bioproject PRJNA496378.
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
# Load metadata
metadat <- read.csv("Raw-data/2021_04_09-data.csv")
metadat2 <- metadat[which(metadat$First.Author.Last.Name=="Gao"),]

##############################################
# This study used the following sequencing approach: 
metadat2$Sequenced.taxa
metadat2$Bacterial.forward.primer
metadat2$Bacterial.reverse.primers
metadat2$Sequencing.mode
metadat2$Sequencing.platform
metadat2$Sequence.processing.approach
# [1] "bacteria"
# [1] "319F"
# [1] "806R"
# [1] "2x300 paired end"
# [1] "Illumina MiSeq"
# [1] "The achieved paired-end reads were assigned to samples based on their unique barcode and truncated by cutting off the barcode and primer sequence. Paired-end reads were merged using FLASH (1.2.8) [23]. Quality filtering on the raw reads was performed by filtering 1) reads containing undetermined nucleotides ratio greater than 5%; 2) reads having the base number of Q (Sequencing quality) ≤10 more than 20%. The specific information of the read number before and after quality filtering is listed in Table 2. Sequences with ≥97% similarity were assigned to the same operational taxonomic units (OTUs) by Usearch."









##############################################
# Directory containing the fastq files after unzipping.
path <- "Raw-data/study-sequences/Gao/seqs-deinterleaved"

##############################################
# Read in the names of the fastq files, and perform some string manipulation 
# to get matched lists of the forward and reverse fastq files.
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = TRUE)) # forward reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



##############################################
# Sequence metadata
seqmet <- read.csv("Raw-data/study-sequences/Gao/SraRunTable.csv")
# select samples taken from monocultures or crop alleys
seqmet2 <- seqmet[which(seqmet$Sample.Name  %in% c("14-yearAgroforestry-cropalley plot1", "14-yearAgroforestry-cropalley plot2", "14-yearAgroforestry-cropalley plot3",
                                                   "5-yearAgroforestry-cropalley plot1", "5-yearAgroforestry-cropalley plot2", "5-yearAgroforestry-cropalley plot3", 
                                                   "9-yearAgroforestry-cropalley plot1", "9-yearAgroforestry-cropalley plot2", "9-yearAgroforestry-cropalley plot3",
                                                   "Cropmonoculture  Plot1", "Cropmonoculture  Plot2", "Cropmonoculture  Plot3")),]
table(seqmet2$cur_land_use) # replication
write.csv(seqmet2, "Raw-data/study-sequences/Gao/SraRunTable_updated.csv")

# 12 total samples
# reduce sample names to those that are in the sequence metadata
keep <- which(sample.names %in% seqmet2$Experiment) # samples to keep
sample.names <- sample.names[keep]
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]


##############################################
# Inspect read quality profiles
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs[c(1:12)]) +# forward reads
  geom_hline(yintercept=30) +
  geom_vline(xintercept=200)
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs[c(1:12)]) + # forward reads
  geom_hline(yintercept=30) +
  geom_vline(xintercept=90)
dev.off()


##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# primers used in the study
fprim <- "ACTCCTACGGGAGGCAGCAG"
rprim <- "GGACTACHVGGGTWTCTAAT"

out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,90), 
                            trimLeft = c(nchar(fprim), nchar(rprim)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) 
# truncLen parameters correspond to quality scores of at least 30
head(out)


##############################################
# Learn the Error Rates

errF <- dada2::learnErrors(filtFs, multithread=TRUE)

errR <- dada2::learnErrors(filtRs, multithread=TRUE)

# It is always worthwhile, as a sanity check if nothing else, to visualize the 
# estimated error rates:

errorplot <- dada2::plotErrors(errF, nominalQ=TRUE)
pdf(paste0(path,"/plots/errorplot_F.pdf"), height=8, width=8)
plot(errorplot)
dev.off()

errorplot <- dada2::plotErrors(errR, nominalQ=TRUE)
pdf(paste0(path,"/plots/errorplot_R.pdf"), height=8, width=8)
plot(errorplot)
dev.off()


##############################################
# Dereplicate

derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names

derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names


##############################################
# Sample Inference

dadaFs <- dada2::dada(derep_forward, err=errF, multithread=TRUE)

dadaRs <- dada2::dada(derep_reverse, err=errR, multithread=TRUE)

# Inspecting the returned dada-class object:
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 352 sequence variants were inferred from 16151 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16



##############################################
# Merge paired reads

mergers <- dada2::mergePairs(dadaFs, derep_forward, 
                             dadaRs, derep_reverse, 
                             verbose=FALSE, justConcatenate = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


##############################################
# Construct sequence table

seqtab <- dada2::makeSequenceTable(mergers)

# row number is your sample number, column number is your ASV number
dim(seqtab)
## [1]   12 54657



# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.

# 270 (this is the sequence length)
# 54657 (this is the number of sequences)




##############################################
# Remove chimeras

seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
## Identified 49419 bimeras out of 54657 input sequences.

dim(seqtab.nochim)
## [1]   12 5238

dim(seqtab.nochim)[2]/dim(seqtab)[2]
## 0.09583402: proportion of sequence variants that are not chimeric
sum(seqtab.nochim)/sum(seqtab)
## 0.6444856: proportion of sequences that are not chimeric






##############################################
# Track reads through the pipeline

getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               round(rowSums(seqtab.nochim)/out[,1]*100, 1))
colnames(track) <- c("input", 
                     "filtered", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", 
                     "nonchim",
                     "final_perc_reads_retained")
rownames(track) <- sample.names
head(track)
track <- as.data.frame(track)
write.csv(track, paste0(path,"/track.csv"))

# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)

track <- read.csv(paste0(path,"/track.csv"))
sum(track$nonchim) # total number of reads: 173675
173675/dim(track)[1]



##############################################
# Save the files
seqtab.nochim <- t(seqtab.nochim)
write.csv(seqtab.nochim, file=paste0(path,"/seqtabs/seqtab1.csv"))
#save(list=ls(), file=list.csv") 





##############################################
# Summary tables
seqtab.nochim <- t(seqtab.nochim)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(path,"/seqtabs/ASVs1.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(path,"/seqtabs/ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)








### read in output files from DADA2 sequence processing code
### these files have the asv sequences in the first column and the
### sequence abundances for each sample in the columns thereafter
# bacteria
bac <- read.csv(paste0(path,"/seqtabs/seqtab1.csv"), row.names=1)


##### determine if any asvs were detected in pcr blank
# no pcr blank in this dataset
#unique(bac$PCRblank.16S); which(bac$PCRblank.16S>0) 


##### remove pcr blank
#bac <- bac[,-81]



#### remove asvs with less than 10 reads across samples

library(readr) 
seqtab.nochim <- read_csv(paste0(path,"/seqtabs/seqtab1.csv"))
bac <- seqtab.nochim
bac <- as.data.frame(bac)
rownames(bac) <- seqtab.nochim$X1
bac <- bac[,-1]
#bac <- t(bac)
bac <- as.data.frame(bac)


# bacteria
bac$sum <- rowSums(bac)
length(bac$sum) # there are 4936 asvs total
length(which(bac$sum<10)) # there are 2008 asvs with <10 reads
length(which(bac$sum>=10)) # there are 2928 asvs with 10 or greater reads
length(which(bac$sum>=10))/length(bac$sum) # proportion of asvs retained: 0.5931929
bac2 <- bac[-which(bac$sum<10),]
hist(bac2$sum); range(bac2$sum) # check
write.csv(bac2, paste0(path,"/seqtabs/seqtab2_10readsormore.csv"))
bacfasta <- ape::read.FASTA(paste0(path,"/seqtabs/ASVs1.fa"))
bacfasta2 <- bacfasta[-which(bac$sum<10)]
ape::write.FASTA(bacfasta2, paste0(path,"/seqtabs/ASVs2_10readsormore.fa"))









