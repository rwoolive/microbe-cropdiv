


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
# Load metadata
metadat <- read.csv("Raw-data/2021_04_09-data.csv")
metadat2 <- metadat[which(metadat$First.Author.Last.Name=="Ashworth"),]

##############################################
# This study used the following sequencing approach: 
metadat2$Sequenced.taxa
metadat2$Bacterial.forward.primer
metadat2$Bacterial.reverse.primers
metadat2$Fungal.forward.primer
metadat2$Fungal.reverse.primers
metadat2$Sequencing.mode
metadat2$Sequencing.platform
metadat2$Sequence.processing.approach





##############################################
# Directory containing the fastq files after unzipping.
# bacteria
path <- "Raw-data/study-sequences/Ashworth/ash fastq files"
subset <- "/recm-year-2-cccc"

##############################################
# Read in the names of the fastq files, and perform some string manipulation 
# to get matched lists of the forward and reverse fastq files.

# bacteria
fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE)) # reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



##############################################
# Sequence metadata
seqmet <- read.csv("Raw-data/study-sequences/Ashworth/Ash-Metadata.csv")
# select RECM, year 2, corn>corn>corn>corn sequence
seqmet2 <- seqmet[which(seqmet$Location == "RECM" & seqmet$Yr=="2" & seqmet$SequenceLet=="C>C>C>C"),]
reps <- table(seqmet2$SequenceLet, seqmet2$Cover, seqmet2$Yr) # replication
reps
write.csv(reps, paste0(path, subset, "/replication.csv"))


# reduce sample names to those that are in the sequence metadata
keep <- which(sample.names %in% seqmet2$Sequence) # samples to keep
sample.names <- sample.names[keep]
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]


##############################################
# Inspect read quality profiles

# set quality score to 30
qs <- 30

# bacteria
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path, subset, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs[c(1:12)]) + # forward reads
  geom_hline(yintercept=qs) + # quality score
  geom_vline(xintercept=210) # adjust according to your sequences
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path, subset, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs[c(1:12)]) + # reverse reads
  geom_hline(yintercept=qs) +
  geom_vline(xintercept=160) # adjust according to your sequences
dev.off()

# save these values for when we trim
trim <- c(210, 160)



##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# bacteria:
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, subset, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, subset, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# primers used in the study
fprim <- "CCTACGGGNGGCWGCAG"
rprim <- "GACTACHVGGGTATCTAATCC"
# these primers amplify the v3-v4 16S gene region, which has an amplicon size of 
# 464 bp

out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trim, # truncLen parameters correspond to quality scores of at least 30
                            trimLeft = c(nchar(fprim), nchar(rprim)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) 
head(out)






##############################################
# Learn the Error Rates

# bacteria
errF <- dada2::learnErrors(filtFs, multithread=TRUE)
errR <- dada2::learnErrors(filtRs, multithread=TRUE)
# It is always worthwhile, as a sanity check if nothing else, to visualize the 
# estimated error rates:
errorplot <- dada2::plotErrors(errF, nominalQ=TRUE)
pdf(paste0(path,subset, "/plots/errorplot_F.pdf"), height=8, width=8)
plot(errorplot)
dev.off()
errorplot <- dada2::plotErrors(errR, nominalQ=TRUE)
pdf(paste0(path,subset, "/plots/errorplot_R.pdf"), height=8, width=8)
plot(errorplot)
dev.off()




##############################################
# Dereplicate

# bacteria
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names




##############################################
# Sample Inference

# bacteria
dadaFs <- dada2::dada(derep_forward, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(derep_reverse, err=errR, multithread=TRUE)
# Inspecting the returned dada-class object:
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 3511 sequence variants were inferred from 58354 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16






##############################################
# Merge paired reads

# bacteria
mergers <- dada2::mergePairs(dadaFs, derep_forward, 
                             dadaRs, derep_reverse, 
                             verbose=FALSE, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



##############################################
# Construct sequence table


# bacteria
seqtab <- dada2::makeSequenceTable(mergers)
# row number is your sample number, column number is your ASV number
dim(seqtab)
## [1]     48 473110
# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.
# 342 (this is the sequence length)
# 473110 (this is the number of asvs)



##############################################
# Remove chimeras

# bacteria
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", 
                                           multithread=TRUE, verbose=TRUE)
## Identified 249287 bimeras out of 282060 input sequences.
dim(seqtab.nochim)
## [1]   30 32773
1-dim(seqtab.nochim)[2]/dim(seqtab)[2]
## 0.8838084: proportion of sequence variants that are not chimeric
1-sum(seqtab.nochim)/sum(seqtab)
## 0.2293771: proportion of sequences that are not chimeric





##############################################
# Track reads through the pipeline


# bacteria
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
write.csv(track, paste0(path,subset, "/track.csv"))
mean(track$final_perc_reads_retained) # 56.94333% reads retained
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)



#### rename to keep track of different subsets
out_R2 <- out 
seqmet2_R2cccc <- seqmet2
filtFs_R2cccc <- filtFs
filtRs_R2cccc <- filtRs
fnFs_R2cccc <- fnFs
fnRs_R2cccc <- fnRs
sample.names_R2cccc <- sample.names
errF_R2cccc <- errF
errR_R2cccc <- errR
derep_forward_R2cccc <- derep_forward
derep_reverse_R2cccc <- derep_reverse
dadaFs_R2cccc <- dadaFs
dadaRs_R2cccc <- dadaRs
mergers_R2cccc <- mergers
seqtab_R2cccc <- seqtab
seqtab.nochim_R2cccc <- seqtab.nochim
subset_R2cccc <- subset

# then remove original files from the environment
rm(out, seqmet2, filtFs, filtRs, fnFs, fnRs, sample.names, errF, errR,
   derep_forward, derep_reverse, dadaFs, dadaRs, mergers, seqtab, 
   seqtab.nochim, subset)



##############################################
# Save the files

# bacteria
seqtab.nochim_R2cccc <- t(seqtab.nochim_R2cccc)
path <- "Raw-data/study-sequences/Ashworth/ash fastq files"
subset <- "/recm-year-2-cccc"
write.csv(seqtab.nochim_R2cccc, file=paste0(path,subset, "/seqtabs/seqtab1.csv"))
save(list=ls(), file=paste0(path,"list-recm-year-2-cccc.RData"))










