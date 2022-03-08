


##############################################
# DADA2: https://benjjneb.github.io/dada2/index.html



##############################################
# Starting point
# Raw sequence data in the form of compressed .fastq files 
# These files were obtained from NCBI Sequence Read Archive (SRA) database (accession numbers SRP155104)
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
metadat2 <- metadat[which(metadat$First.Author.Last.Name=="Song"),]

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
# [1] "bacteria, archaea, fungi"
# [1] "515F"
# [1] "907R"
# [1] "ITS1F"
# [1] "ITS1R"
# [1] "2×250 paired end cycle"
# [1] "Illumina MiSeq"
# [1] "Fastx was used to remove bases with end masses lower than Q15; the paired sequences obtained by double-end sequencing were merged into a longer sequence by FLASH in an overlapping relationship; and the obtained splicing sequence (merge) was obtained. Cutadapt software was used to delete themiscellaneous sequence at the end of the short sequence, and the resulting sequence was filtered with USEARCH software to obtain the valid sequence (Primer Filter). The singleton sequence was removed by USEARCH software, and the reference standard was compared to remove the chimeric sequence. The UPARSE clusteringmethod was performed tomerge the quality-controlled sequences (97-degree similarity between sequences) into an operational taxonomic unit (OTU) cell for subsequent OTU classification annotations (Edgar, 2013)."






##############################################
# Directory containing the fastq files after unzipping.
# bacteria
path <- "Raw-data/study-sequences/Song/seqs/bac"
# fungi
path2 <- "Raw-data/study-sequences/Song/seqs/fun"

##############################################
# Read in the names of the fastq files, and perform some string manipulation 
# to get matched lists of the forward and reverse fastq files.

# bacteria
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = TRUE)) # reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# fungi
fnFs2 <- sort(list.files(path2, pattern="R1.fastq.gz", full.names = TRUE)) # forward reads
fnRs2 <- sort(list.files(path2, pattern="R2.fastq.gz", full.names = TRUE)) # reverse reads
sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)


##############################################
# Inspect read quality profiles

# set quality score to 30
qs <- 30

# bacteria
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs[c(1:12)]) + # forward reads
  geom_hline(yintercept=qs) + # quality score
  geom_vline(xintercept=240) # adjust according to your sequences
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs[c(1:12)]) + # reverse reads
  geom_hline(yintercept=qs) +
  geom_vline(xintercept=150)
dev.off()

# save these values for when we trim
trim <- c(240, 150)

# fungi
# We start by visualizing the quality profiles of the forward reads:
pdf(paste0(path2, "/plots/qscores_F.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnFs2[c(1:12)]) + # forward reads
  geom_hline(yintercept=qs) + # quality score
  geom_vline(xintercept=225) # adjust according to your sequences
dev.off()

# Now we visualize the quality profile of the reverse reads:
pdf(paste0(path2, "/plots/qscores_R.pdf"), height=8, width=8)
dada2::plotQualityProfile(fnRs2[c(1:12)]) + # reverse reads
  geom_hline(yintercept=qs) +
  geom_vline(xintercept=200)
dev.off()

# save these values for when we trim
trim2 <- c(225, 200)


##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# bacteria:
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# primers used in the study
fprim <- "GTGCCAGCMGCCGCGG"
rprim <- "CCGTCAATTCMTTTRAGTTT"
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trim, # truncLen parameters correspond to quality scores of at least 30
                            trimLeft = c(nchar(fprim), nchar(rprim)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) 
head(out)

# fungi:
# Place filtered files in filtered/ subdirectory
filtFs2 <- file.path(path2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path2, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2
# primers used in the study
fprim2 <- "CTTGGTCATTTAGAGGAAGTAA"
rprim2 <- "GCTGCGTTCTTCATCGATGC"
out2 <- dada2::filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=trim2, # truncLen parameters correspond to quality scores of at least 30
                            trimLeft = c(nchar(fprim2), nchar(rprim2)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) 
head(out2)




##############################################
# Learn the Error Rates

# bacteria
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

# fungi
errF2 <- dada2::learnErrors(filtFs2, multithread=TRUE)
errR2 <- dada2::learnErrors(filtRs2, multithread=TRUE)
# It is always worthwhile, as a sanity check if nothing else, to visualize the 
# estimated error rates:
errorplot2 <- dada2::plotErrors(errF2, nominalQ=TRUE)
pdf(paste0(path2,"/plots/errorplot_F.pdf"), height=8, width=8)
plot(errorplot2)
dev.off()
errorplot2 <- dada2::plotErrors(errR2, nominalQ=TRUE)
pdf(paste0(path2,"/plots/errorplot_R.pdf"), height=8, width=8)
plot(errorplot2)
dev.off()



##############################################
# Dereplicate

# bacteria
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# fungi
derep_forward2 <- derepFastq(filtFs2, verbose=TRUE)
names(derep_forward2) <- sample.names2
derep_reverse2 <- derepFastq(filtRs2, verbose=TRUE)
names(derep_reverse2) <- sample.names2


##############################################
# Sample Inference

# bacteria
dadaFs <- dada2::dada(derep_forward, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(derep_reverse, err=errR, multithread=TRUE)
# Inspecting the returned dada-class object:
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 4502 sequence variants were inferred from 129841 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


# fungi
dadaFs2 <- dada2::dada(derep_forward2, err=errF2, multithread=TRUE)
dadaRs2 <- dada2::dada(derep_reverse2, err=errR2, multithread=TRUE)
# Inspecting the returned dada-class object:
dadaFs2[[1]]
# dada-class: object describing DADA2 denoising results
# 731 sequence variants were inferred from 48987 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16





##############################################
# Merge paired reads

# bacteria
mergers <- dada2::mergePairs(dadaFs, derep_forward, 
                             dadaRs, derep_reverse, 
                             verbose=FALSE, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# fungi
mergers2 <- dada2::mergePairs(dadaFs2, derep_forward2, 
                             dadaRs2, derep_reverse2, 
                             verbose=FALSE, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers2[[1]])


##############################################
# Construct sequence table


# bacteria
seqtab <- dada2::makeSequenceTable(mergers)
# row number is your sample number, column number is your ASV number
dim(seqtab)
## [1]     12 135292
# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.
# 364 (this is the sequence length)
# 135292 (this is the number of asvs)


# fungi
seqtab2 <- dada2::makeSequenceTable(mergers2)
# row number is your sample number, column number is your ASV number
dim(seqtab2)
## [1]    12 7950
# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab2)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.
# 393 (this is the sequence length)
# 7950 (this is the number of asvs)




##############################################
# Remove chimeras

# bacteria
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
## Identified 105162 bimeras out of 135292 input sequences.
dim(seqtab.nochim)
## [1]   12 30130
1-dim(seqtab.nochim)[2]/dim(seqtab)[2]
## 0.7772965: proportion of sequence variants that are not chimeric
1-sum(seqtab.nochim)/sum(seqtab)
## 0.1094492: proportion of sequences that are not chimeric

# fungi
seqtab.nochim2 <- dada2::removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
## Identified 1170 bimeras out of 7950 input sequences.
dim(seqtab.nochim2)
## [1]   12 6780
1-dim(seqtab.nochim2)[2]/dim(seqtab2)[2]
## 0.1471698: proportion of sequence variants that are not chimeric
1-sum(seqtab.nochim2)/sum(seqtab2)
## 0.004665494: proportion of sequences that are not chimeric






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
write.csv(track, paste0(path,"/track.csv"))
mean(track$final_perc_reads_retained) # 75.40833% reads retained
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)


track <- read.csv(paste0(path,"/track.csv"))
sum(track$nonchim) # total number of reads: 2339233
2339233/dim(track)[1]


# fungi
getN <- function(x) sum(dada2::getUniques(x))
track2 <- cbind(out2, 
               sapply(dadaFs2, getN), 
               sapply(dadaRs2, getN), 
               sapply(mergers2, getN), 
               rowSums(seqtab.nochim2),
               round(rowSums(seqtab.nochim2)/out2[,1]*100, 1))
colnames(track2) <- c("input", 
                     "filtered", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", 
                     "nonchim",
                     "final_perc_reads_retained")
rownames(track2) <- sample.names
head(track2)
track2 <- as.data.frame(track2)
write.csv(track2, paste0(path2,"/track.csv"))
mean(track2$final_perc_reads_retained) # .75% reads retained
# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)

track <- read.csv(paste0(path2,"/track.csv"))
sum(track$nonchim) # total number of reads: 2894805
2894805/dim(track)[1]


##############################################
# Save the files

# bacteria
seqtab.nochim <- t(seqtab.nochim)
write.csv(seqtab.nochim, file=paste0(path,"/seqtabs/seqtab1.csv"))
#save(list=ls(), file=list.csv") 

# fungi
seqtab.nochim2 <- t(seqtab.nochim2)
write.csv(seqtab.nochim2, file=paste0(path2,"/seqtabs/seqtab1.csv"))
#save(list=ls(), file=list.csv") 




##############################################
# Summary tables

# bacteria
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



# fungi
seqtab.nochim2 <- t(seqtab.nochim2)
asv_seqs2 <- colnames(seqtab.nochim2)
asv_headers2 <- vector(dim(seqtab.nochim2)[2], mode="character")
for (i in 1:dim(seqtab.nochim2)[2]) {
  asv_headers2[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta2 <- c(rbind(asv_headers2, asv_seqs2))
write(asv_fasta2, paste0(path2,"/seqtabs/ASVs1.fa"))
# count table:
asv_tab2 <- t(seqtab.nochim2)
row.names(asv_tab2) <- sub(">", "", asv_headers2)
write.table(asv_tab2, paste0(path2,"/seqtabs/ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)



