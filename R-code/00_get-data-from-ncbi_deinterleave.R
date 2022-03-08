# At this point we have identified PRJ numbers for 
# studies that we want to get sequence data from.
# First, go to:
# https://www.ncbi.nlm.nih.gov/Traces/study/?go=home
# Enter the PRJ number (e.g., PRJNA496378 for Gao 2019) 
# into the search box and hit enter.
# This will take you to a page that shows the sequence
# files within that project. 
# Select all runs at the bottom of the page, then download
# the metadata and accession list. Put those files into
# the Raw-data/study-sequences/(Author-Year) directory.
# 
# Now go to https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name
# Copy one experiment number at a time from the SraRunTable 
# (metadata table), paste it into the search box and hit 
# enter. This will take you to a new page. Make sure that the
# "FASTQ" download format is selected, and hit download. 
# Put this file into the Raw-data/study-sequences/(Author-Year)/seqs 
# directory and replace sra_data with the Sample Name from the SRA run 
# selector page. 

# Once you have saved all of the fastq files, you may need to deinterleave
# the paired reads. DADA can only work with fastq.gz files that are separated 
# by forward/reverse read format.
# To deinterleave, install the GGMap software (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
# and put the GGmap tar.gz file onto your desktop
# Go to your terminal and change your working directory to your Desktop
cd Desktop
# Extract info from the BBMap tar.gz
tar -xvzf BBMap_38.90.tar.gz
# Drag all fastq.gz files onto your desktop  
# Split files (one by one) into forward and reverse reads, e.g.:
bbmap/reformat.sh in=SRX4516068.fastq.gz out1=SRX4516068_R1.fastq.gz out2=SRX4516068_R2.fastq.gz







# Now you are ready to process the sequences using DADA2 
# (next R code)



