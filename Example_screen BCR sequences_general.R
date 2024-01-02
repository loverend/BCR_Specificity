#### Stage 1 set up files
## Combine sequences with reference ready for python pipeline 
## Rachael Bashford-Rogers and Lauren Overend 
## October 2022

#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
library(seqinr)
#...................................................................
setwd('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining')
## Outputdir from BCR TCR processing pipeline
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
functionality <- "PRODUCTIVE"
runname <- "Screen_Sepsis"

## PreProcessing: Lets get the CDR3s for all of my samples from IMGT!!!
source("GetCDR3s.R")
get_cdr3s(outputdir, functionality)
### Lets get all the rds for CDR3s
cdr3dir <- paste0(outputdir, "/CDR3_AA")
## CDR3s are agregated within the same individual - not across individuals/timepoints
all_rds <- list.files(cdr3dir, full.name=TRUE)
all_rds <- grep(".rds", all_rds, value=TRUE)
all_rds <- grep("ALL_SAMPLES", all_rds, value=TRUE, invert=TRUE)
all_rds1 <- all_rds
####
new_dir <- paste0(cdr3dir, "/FOR_MINING")
new_dirx <-new_dir
if (!dir.exists(new_dir)) {dir.create(new_dir)}

#...................................................................
# Load Panel of reference sequences 
# If you have not got this file you must first run BuildSequenceReference.R
reference_seqs <- readRDS("HUMAN_PATHOGEN_UK_REFERENCE.rds")
#.....................................................................


##......................................................................................
####### Stage 3 Merge sequences together and output per sample for python analysis  ##
for(i in 1:length(all_rds)){
	sampleid <- basename(all_rds)[i]
	sampleid <- gsub(".rds", "_screen_pathogen_sequences.txt", sampleid)
	screen_seqs <- readRDS(all_rds[i])
	####### Stage 3 filter for long CDR3s
	screen_seqs<-screen_seqs[nchar(screen_seqs)>=5]
	reference_seqs<-reference_seqs[nchar(reference_seqs)>=5]
	screen_seqs = screen_seqs[grep("*", screen_seqs, invert = T, fixed = T)]
	####### write out for python script
	x = rbind(cbind(screen_seqs, "screen"),cbind(reference_seqs, "reference"))
	s <- c(screen_seqs, reference_seqs)
	s <- as.list(s)
	## Save to the file for python analysis 
	output_file = paste0(new_dir, "/", sampleid)
	sampleid <- gsub(".txt", "", sampleid)
	output_fileq = paste0(new_dir, "/", sampleid, "_WITHSAMPLEID.txt")
	output_filef = paste0(new_dir, "/", sampleid, "_screen.fasta")
	write.table(x, file = output_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
	write.fasta(s, names(s), output_filef, open = "w", nbchar = 60, as.string = FALSE)
	write.table(x, file = output_fileq, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".",col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

## Now we make a sample file listing the links to all the files for analysis
all_rds <- list.files(new_dir, full.name=TRUE)
all_rds <- grep("WITHSAMPLEID", all_rds, value=TRUE, invert=TRUE)
all_rds <- grep("fasta", all_rds, value=TRUE, invert=TRUE)

## Setting the sample output names
all_rds <- data.frame(all_rds)
all_rds[,2] <- all_rds[,1]
all_rds[,2] <- gsub("FOR_MINING", "RESULTS", all_rds[,2])
all_rds[,2] <- gsub("screen_pathogen_sequences", "matched_cdr3s", all_rds[,2])

## Make the outputdirectory
new_dir <- gsub("FOR_MINING", "RESULTS", new_dir)
if (!dir.exists(new_dir)) {dir.create(new_dir)}

##saving
outputfile <- paste0(runname, "_sample_spreadsheet_for_mining.txt")
write.table(all_rds, outputfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

###################### Extract in fasta format to run cdhit 
all_rds <- list.files(new_dirx, full.name=TRUE)
all_rds <- grep(".fasta", all_rds, value=TRUE)

## Setting the sample output names
all_rds <- data.frame(all_rds)
all_rds[,2] <- all_rds[,1]
all_rds[,2] <- gsub("FOR_MINING", "CDHIT", all_rds[,2])
all_rds[,2] <- gsub("screen_pathogen_sequences", "cdhit_out", all_rds[,2])
all_rds[,2] <- gsub(".fasta", " ", all_rds[,2])

## Make the outputdirectory
new_dir <- gsub("FOR_MINING", "CDHIT", new_dir)
if (!dir.exists(new_dir)) {dir.create(new_dir)}

##saving
outputfile <- paste0(runname, "_sample_spreadsheet_for_cdhit.txt")
write.table(all_rds, outputfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

##...................................................................................... 
##......................................................................................
## Actually maybe it makes sense to merge files from the same individual!
## Lets make one fasta file per sample - rather than per timepoint 
all_rds <- all_rds1
samples_all <- basename(all_rds)
samples_all <- gsub("_PRODUCTIVE_CDR3AA.rds", "", samples_all)
sample_use <- c()
for(i in 1:length(samples_all)){
	samplex <- samples_all[i]
	if(samplex %like% "HV"){
		a <- str_split_fixed(samplex, "_", 3)
		a <- paste0(a[1], "_", a[,2])
	} else {
		a <- str_split_fixed(samplex, "_", 3)[,1]
	}
	sample_use <- c(sample_use, a)
}
sample_use <- unique(sample_use)

for(i in 1:length(sample_use)){
	sampleid <- sample_use[i]
	files_all <- grep(sampleid, all_rds1, value=TRUE)
	seqs <- c()
	for(x in 1:length(files_all)){
		d <- readRDS(files_all[x])
		seqs <- c(seqs, d)
	}
	screen_seqs <- seqs
	####### Stage 3 filter for long CDR3s
	screen_seqs<-screen_seqs[nchar(screen_seqs)>=5]
	reference_seqs<-reference_seqs[nchar(reference_seqs)>=5]
	screen_seqs = screen_seqs[grep("*", screen_seqs, invert = T, fixed = T)]
	####### write out for python script
	x = rbind(cbind(screen_seqs, "screen"),cbind(reference_seqs, "reference"))
	s <- c(screen_seqs, reference_seqs)
	s <- as.list(s)
	## Save to the file for python analysis 
	output_file = paste0(new_dirx, "/", sampleid, ".txt")
	sampleid <- gsub(".txt", "", sampleid)
	output_fileq = paste0(new_dirx, "/", sampleid, "_ALLSAMPLES_WITHSAMPLEID.txt")
	output_filef = paste0(new_dirx, "/", sampleid, "_ALLSAMPLES_screen.fasta")
	write.table(x, file = output_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
	write.fasta(s, names(s), output_filef, open = "w", nbchar = 60, as.string = FALSE)
	write.table(x, file = output_fileq, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".",col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

## Now we make a sample file listing the links to all the files for analysis
all_rds <- list.files(new_dirx, full.name=TRUE)
all_rds <- grep("WITHSAMPLEID", all_rds, value=TRUE, invert=TRUE)
all_rds <- grep("_ALLSAMPLES", all_rds, value=TRUE)

## Setting the sample output names
all_rds <- data.frame(all_rds)
all_rds[,2] <- all_rds[,1]
all_rds[,2] <- gsub("FOR_MINING", "RESULTS_PERINDIVIDUAL", all_rds[,2])
all_rds[,2] <- gsub("screen_pathogen_sequences", "matched_cdr3s", all_rds[,2])

## Make the outputdirectory
new_dir <- gsub("FOR_MINING", "RESULTS_PERINDIVIDUAL", new_dirx)
if (!dir.exists(new_dir)) {dir.create(new_dir)}

##saving
outputfile <- paste0(runname, "_sample_spreadsheet_for_mining_PERINDIVIDUAL.txt")
write.table(all_rds, outputfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

###################### Extract in fasta format to run cdhit 
all_rds <- list.files(new_dirx, full.name=TRUE)
all_rds <- grep(".fasta", all_rds, value=TRUE)
all_rds <- grep("_ALLSAMPLES", all_rds, value=TRUE)

## Setting the sample output names
all_rds <- data.frame(all_rds)
all_rds[,2] <- all_rds[,1]
all_rds[,2] <- gsub("FOR_MINING", "CDHIT_PERINDIVIDUAL", all_rds[,2])
all_rds[,2] <- gsub("screen_pathogen_sequences", "cdhit_out", all_rds[,2])
all_rds[,2] <- gsub(".fasta", " ", all_rds[,2])

## Make the outputdirectory
new_dir <- gsub("FOR_MINING", "CDHIT_PERINDIVIDUAL", new_dirx)
if (!dir.exists(new_dir)) {dir.create(new_dir)}

##saving
outputfile <- paste0(runname, "_sample_spreadsheet_for_cdhit_PERINDIVIDUAL.txt")
write.table(all_rds, outputfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
##......................................................................................