### Function to get CDR3s from IMGT output and save as a file per individual 
## Lauren Overend
#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
#functionality <- "PRODUCTIVE"
 
library(stringr)
 
 get_cdr3s <- function(outputdir, functionality){
	directory_imgt <- paste0(outputdir, "/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
	files <- list.files(directory_imgt, full.name=TRUE)
	if(functionality=="PRODUCTIVE"){
		files <- grep("_productive_5_AA-sequences.txt", files, value=TRUE)
	}
	
	### Aggregating cdr3s within the same file (timepoint) not across timepoints 
	all_cdr3s <- c()
	for(i in 1:length(files)){
		sampleID <- basename(files[i])
		seqs <- read.delim(files[i], header=FALSE, sep="\t")
		cdr3s <- seqs[, c(1,2,15)]
		cdr3s$seqid <- str_split_fixed(cdr3s$V2, "__", 2)[,1]
		cdr3s <- cdr3s[, c("seqid", "V15")]
		## want to aggregate into unique cdr3s
		cdr3sx <- aggregate(data=cdr3s,seqid~V15,FUN=paste,collapse=",")
		cdr3s <- cdr3sx
		cdr3s$cdr3count <- (str_count(cdr3s$seqid, ",")+1)
		cdr3s$index <- rownames(cdr3s)
		
		###########################################
		sampleID <- gsub("IMGT_", "", sampleID)
		sampleID <- gsub("_productive_5_AA-sequences.txt", "", sampleID)
		##########################################
		cdr3s$index	<- paste0(sampleID, ":", cdr3s$index)
		cdr3dir <- paste0(outputdir, "/CDR3_AA")
		if (!dir.exists(cdr3dir)) {dir.create(cdr3dir)}
		filename <- paste0(cdr3dir, "/", sampleID, "_", functionality, "_CDR3AA.txt") 
		write.table(cdr3s, filename, sep="\t", row.names=FALSE)	
		#############################
		cdr3 <- cdr3s[c(1,4)]
		rownames(cdr3) <- cdr3[,2]
		cdr3[,2] <- NULL
		cdr3 <- cdr3[,1]
		names(cdr3) <- cdr3s$index
		all_cdr3s <- c(all_cdr3s, cdr3)
		filename <- paste0(cdr3dir, "/", sampleID, "_", functionality, "_CDR3AA.rds") 
		saveRDS(cdr3, filename)
	}
	filename <- paste0(cdr3dir, "/ALL_SAMPLES_", functionality, "_CDR3AA.txt") 
	write.table(all_cdr3s, filename, sep="\t", col.names=FALSE)
	filename <- paste0(cdr3dir, "/ALL_SAMPLES_", functionality, "_CDR3AA.RDS") 
	saveRDS(all_cdr3s, filename)
}