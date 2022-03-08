


##############################################
# DADA2: https://benjjneb.github.io/dada2/index.html



##############################################
# Starting point
# Raw sequence data in the form of compressed .fastq files 
# These files were obtained from National Center for Biotechnology Information under Bioproject ID: PRJNA481943
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
metadat2 <- metadat[which(metadat$First.Author.Last.Name=="Cloutier" & metadat$Publication.Year==2020),]

##############################################
# This study used the following sequencing approach: 
metadat2$Sequenced.taxa
metadat2$Fungal.forward.primer
metadat2$Fungal.reverse.primers
metadat2$Sequencing.mode
metadat2$Sequencing.platform
metadat2$Sequence.processing.approach





##############################################
# Directory containing the fastq files after unzipping.
# fungi
path <- "Raw-data/study-sequences/Cloutier/seqs-deinterleaved"

##############################################
# Read in the names of the fastq files, and perform some string manipulation 
# to get matched lists of the forward and reverse fastq files.

# fungi
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE)) # reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)



##############################################
# Sequence metadata
seqmet <- read.csv("Raw-data/study-sequences/Cloutier/SraRunTable.csv")

# make a factor for cover type
seqmet$Cover <- as.factor(sapply(strsplit(basename(seqmet$Library.Name), " CC Block"), `[`, 1))
seqmet$Cover <- factor(seqmet$Cover, levels(seqmet$Cover)[c(5,3,4,6,7,8,9,1,2)])

# make a factor for replicate
seqmet$Rep <- sapply(strsplit(basename(seqmet$Library.Name), " CC Block "), `[`, 2)
seqmet$Rep <- as.factor(sapply(strsplit(basename(seqmet$Rep), " "), `[`, 1))

# make a factor for season
seqmet$Season <- sapply(strsplit(basename(seqmet$Library.Name), " CC Block "), `[`, 2)
seqmet$Season <- as.factor(sapply(strsplit(basename(seqmet$Season), " "), `[`, 2))

table(seqmet$Cover, seqmet$Season) # replication

# list of sequence IDs to keep
length(seqmet$Experiment %in% sample.names)

# 72 total samples (9 covers x 2 seasons x 4 reps)



##############################################
# Inspect read quality profiles

# set quality score to 30
qs <- 30

# fungi
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs[c(1:12)]) + # forward reads
  geom_hline(yintercept=qs) + # quality score
  geom_vline(xintercept=250) # adjust according to your sequences
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs[c(1:12)]) + # reverse reads
  geom_hline(yintercept=qs) +
  geom_vline(xintercept=250) # adjust according to your sequences
dev.off()

# save these values for when we trim
trim <- c(250, 250)



##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# fungi:
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# primers used in the study
fprim <- "CTTGGTCATTTAGAGGAAGTAA"
rprim <- "GCTGCGTTCTTCATCGATGC"
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trim, # truncLen parameters correspond to quality scores of at least 30
                            trimLeft = c(nchar(fprim), nchar(rprim)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) 
head(out)






##############################################
# Learn the Error Rates

# fungi
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

# fungi
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names



##############################################
# Sample Inference

# fungi
dadaFs <- dada2::dada(derep_forward, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(derep_reverse, err=errR, multithread=TRUE)
# Inspecting the returned dada-class object:
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 477 sequence variants were inferred from 14042 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16






##############################################
# Merge paired reads

# fungi
mergers <- dada2::mergePairs(dadaFs, derep_forward, 
                             dadaRs, derep_reverse, 
                             verbose=FALSE, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



##############################################
# Construct sequence table


# fungi
seqtab <- dada2::makeSequenceTable(mergers)
# row number is your sample number, column number is your ASV number
dim(seqtab)
## [1]     72 392080
# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.
# 468 (this is the sequence length)
# 392080 (this is the number of asvs)



##############################################
# Remove chimeras

# fungi
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
## Identified 283859 bimeras out of 392080 input sequences.
dim(seqtab.nochim)
## [1]   72 108221
1-dim(seqtab.nochim)[2]/dim(seqtab)[2]
## 0.7239824: proportion of sequence variants that are not chimeric
1-sum(seqtab.nochim)/sum(seqtab)
## 0.4915207: proportion of sequences that are not chimeric





##############################################
# Track reads through the pipeline


# fungi
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
mean(track$final_perc_reads_retained) # 41.66111% reads retained
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)


track <- read.csv(paste0(path,"/track.csv"))
sum(track$nonchim) # total number of reads: 1480374
1480374/dim(track)[1]

##############################################
# Save the files

# fungi
seqtab.nochim <- t(seqtab.nochim)
write.csv(seqtab.nochim, file=paste0(path,"/seqtabs/seqtab1.csv"))
#save(list=ls(), file=list.csv") 



##############################################
# Summary tables

# fungi
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
# fungi
fun <- read.csv(paste0(path,"/seqtabs/seqtab1.csv"), row.names=1)


##### determine if any asvs were detected in pcr blank
# no pcr blank in this dataset
#unique(fun$PCRblank.16S); which(fun$PCRblank.16S>0) 


##### remove pcr blank
#fun <- fun[,-81]



#### remove asvs with less than 10 reads across samples



# fungi
fun$sum <- rowSums(fun)
length(fun$sum) # there are 108221 asvs total
length(which(fun$sum<10)) # there are 93520 asvs with <10 reads
length(which(fun$sum>=10)) # there are 14701 asvs with 10 or greater reads
length(which(fun$sum>=10))/length(fun$sum) # proportion of asvs retained: 0.530103
fun2 <- fun[-which(fun$sum<10),]
hist(fun2$sum); range(fun2$sum) # check
write.csv(fun2, paste0(path,"/seqtabs/seqtab2_10readsormore.csv"))
funfasta <- ape::read.FASTA(paste0(path,"/seqtabs/ASVs1.fa"))
funfasta2 <- funfasta[-which(fun$sum<10)]
ape::write.FASTA(funfasta2, paste0(path,"/seqtabs/ASVs2_10readsormore.fa"))


