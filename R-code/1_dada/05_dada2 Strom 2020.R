


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
metadat2 <- metadat[which(metadat$First.Author.Last.Name=="Strom" & metadat$Publication.Year==2020),]

##############################################
# This study used the following sequencing approach: 
metadat2$Sequenced.taxa
metadat2$fungil.forward.primer
metadat2$fungil.reverse.primers
metadat2$Fungal.forward.primer
metadat2$Fungal.reverse.primers
metadat2$Sequencing.mode
metadat2$Sequencing.platform
metadat2$Sequence.processing.approach





##############################################
# Directory containing the fastq files after unzipping.
# fungi
path <- "Raw-data/study-sequences/Strom/seqs"

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
seqmet <- read.csv("Raw-data/study-sequences/Strom/SraRunTable.csv")
# select Cc, Ss, Ca, and Sa rotation treatments, both years (2015, 2016), all seasons (fall & spring & summer), & exclude technical replicates
seqmet2 <- seqmet[which(seqmet$crop_rotation %in% c("Cc", "Ss", "Ca", "Sa") & is.na(seqmet$replicate)==TRUE),]
table(seqmet2$crop_rotation, seqmet2$block, seqmet2$year) # replication
# list of sequence IDs to keep
seqmet2$Experiment %in% sample.names
keep <- sample.names %in%  seqmet2$Experiment
sample.names <- sample.names[keep]

length(sample.names)
# 96 total samples

fnFs <- fnFs[keep]
fnRs <- fnRs[keep]

##############################################
# Inspect read quality profiles

# set quality score to 30
qs <- 30

# fungi
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs[c(1:12)]) + # forward reads
  geom_hline(yintercept=qs) + # quality score
  geom_vline(xintercept=140) # adjust according to your sequences
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs[c(1:12)]) + # reverse reads
  geom_hline(yintercept=qs) +
  geom_vline(xintercept=120) # adjust according to your sequences
dev.off()

# save these values for when we trim
trim <- c(140, 120)



##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# fungi:
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# primers used in the study (not stated)
fprim <- ""
rprim <- ""
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trim, # truncLen parameters correspond to quality scores of at least 30
                            #trimLeft = c(nchar(fprim), nchar(rprim)),
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
# 382 sequence variants were inferred from 15019 input unique sequences.
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
## [1]     16 17215
# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.
# 270 (this is the sequence length)
# 17215 (this is the number of asvs)



##############################################
# Remove chimeras

# fungi
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
## Identified 11003 bimeras out of 17215 input sequences.
dim(seqtab.nochim)
## [1]   16 6212
1-dim(seqtab.nochim)[2]/dim(seqtab)[2]
## 0.6391519: proportion of sequence variants that are not chimeric
1-sum(seqtab.nochim)/sum(seqtab)
## 0.02990887: proportion of sequences that are not chimeric





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
mean(track$final_perc_reads_retained) # *** 88.9875% reads retained
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)


track <- read.csv(paste0(path,"/track.csv"))
sum(track$nonchim) # total number of reads: 7697025
7697025/dim(track)[1]



##############################################
# Save the files

# fungi
seqtab.nochim <- t(seqtab.nochim)
write.csv(seqtab.nochim, file=paste0(path,"/seqtabs/seqtab1.csv"))
save(list=ls(), file=paste0(path, "_list.RData")) 




##############################################
# Summary tables
load(file=paste0(path, "_list.RData"))

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
length(fun$sum) # there are 18358 asvs total
length(which(fun$sum<10)) # there are 9068 asvs with <10 reads
length(which(fun$sum>=10)) # there are 9290 asvs with 10 or greater reads
length(which(fun$sum>=10))/length(fun$sum) # proportion of asvs retained: 0.5060464
fun2 <- fun[-which(fun$sum<10),]
hist(fun2$sum); range(fun2$sum) # check
write.csv(fun2, paste0(path,"/seqtabs/seqtab2_10readsormore.csv"))
funfasta <- ape::read.FASTA(paste0(path,"/seqtabs/ASVs1.fa"))
funfasta2 <- funfasta[-which(fun$sum<10)]
ape::write.FASTA(funfasta2, paste0(path,"/seqtabs/ASVs2_10readsormore.fa"))


