## Code to format the output of cdhit into a dataframe for easier analysis
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
## November 2022
## Some code adapted from ttps://rpubs.com/rmurdoch/cdhit_to_mapping_file


library(dplyr)
library(tibble)
library(ggplot2)
library(viridisLite)
library(stringr)
library(data.table)
library(phangorn)
library(compbio4all)
library(rentrez)
library(seqinr)
library(msa)
library(ggmsa)
library(phangorn)
library(ape)
library(ggtree)
library(cowplot)
library(foreach)
library(doParallel)

format_cdhit <- function(results, i) {	
	print(i)
	sampleid <- basename(results[i])
	sampleid <- gsub("_ALLSAMPLES_screen.fasta.clstr", "", sampleid)
	sampleid <- gsub("_ALLSAMPLES_screen.clstr", "", sampleid)
	## Check for results
	info = file.info(results[i])
	if(info$size==0){
		next 
	}
	clstr <- read.csv(results[i], sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
	clstr2 <- clstr
	n = nrow(clstr)
	x = 0
	numbers_only <- function(x) !grepl("\\D", x)
	for (row in c(1:n)) {
	  if (numbers_only(clstr2[row,1]) == TRUE) {
		clstr2[row,1] <- x}
	  else {NULL}
	  x <- clstr2[row,1]
	}
	#head(clstr2)
	clstr.sums <- data.frame(dplyr::count(clstr2,V1))
	clstr.sums$n <- clstr.sums$n-1
	#head(clstr.sums)
	switch <- clstr.sums[1,2]
	clstr3 <- cbind(clstr2[1], clstr)
	#clstr3[c((switch-5):(switch+5)),]
	clstr4 <- clstr2[-which(clstr2$V2 == ""), ]
	clstr5 <- clstr4
	clstr5[] <- lapply(clstr5, gsub, pattern='>', replacement='')
	clstr5.2 <- data.frame(str_split_fixed(clstr5$V2, "aa, ", 2))
	clstr5.3 <- data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
	clstr6 <- cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
	colnames(clstr6) <- c("cluster","aa","gene","stat")
	
	## save nice format 
	clstr6$type <- NA
	clstr6$type[clstr6$gene %like% "cdr3"] <- "ref_cdr3"
	clstr6$type[!clstr6$gene %like% "cdr3"] <- sampleid
    clstr6$sample <- sampleid
   	write.table(clstr6, paste0(outputdir, "Hits_", sampleid, ".txt"), sep="\t")
	
	### Still has the cluster up to here 
	## are there any clusters made up of both reference and sample seq 
	f<- table(clstr6$cluster, clstr6$type)
	f <- as.data.frame.matrix(f)
	## Need at least one reference hit and one sample hit!!!
	f <- f[f[,1]>=1,]
	f <- f[f[,2]>=1,]
	
	## list of clusters to keep
	f<- (data.frame(f))
	if(dim(f)[1] > 0){
		print(sampleid)
		print(i)
		cluster_keep <-rownames(f)
		clstr6 <- clstr6[clstr6$cluster %in% cluster_keep,]
		## We want to get the sequence that matches this cluster 
		seqs <- read.delim(paste0(base_dir, "FOR_MINING/", sampleid, "_ALLSAMPLES_WITHSAMPLEID.txt"), header=TRUE)
		seqs <- data.frame(seqs)
		seqs$id <- rownames(seqs)
		seqs$X<- NULL
		clstr6x <- merge(clstr6, seqs, by.x="gene", by.y="id", all.x=TRUE)
		clstr6x <- merge(clstr6x, reference, by.x="gene", by.y="cdr3_id", all.x="TRUE")
		plotdir <- paste0(outputdir, "PLOTS")
		if (!dir.exists(plotdir)) {dir.create(plotdir)}
		pdf(paste0(plotdir, "/CDR3matches_aligned_", sampleid, ".pdf"), width=8, height=4)
		## Align the sequences and plot 
		for(x in 1:length(unique(clstr6x$cluster))){
			new_cluster <- clstr6x[clstr6x$cluster==unique(clstr6x$cluster)[x],]
			pathogen <- unique(new_cluster$pathogen)
			pathogen <- pathogen[!is.na(pathogen)] 
			pathogen<- paste(pathogen, collapse=" ")
			seqs_aa <- new_cluster$screen_seqs
			names(seqs_aa) <- new_cluster$gene 
			aa_ss <- Biostrings::AAStringSet(seqs_aa)
			##d <- as.dist(stringDist(aa_ss, method = "levenshtein")/width(aa_ss)[1])
			dio1s_msa <- msa(aa_ss)
			class(dio1s_msa) <- "AAMultipleAlignment"
			p1 <- ggmsa(dio1s_msa, seq_name = TRUE, color = "Chemistry_AA", show.legend = TRUE) + ggtitle(pathogen)
			#tree <- bionj(d)
			#p2 <- ggtree(tree, ignore.negative.edge=TRUE) + geom_tiplab(size=1)
			#plot(plot_grid(p1, p2, ncol=1))
			plot(p1)
		}
		dev.off()
		return(clstr6x)
	}
}