## Code to subsample BCR reads and assign pathogen hit frequency and counts
# Lauren Overend
# lauren.overend@oriel.ox.ac.uk
# Nov 2022

subsample_cdhits <- function(all_cdhits, i, base_dir, bcr_sequences_iso, sub_depth, plot_dir){
	
	## Setting up Data
	datai <- read.delim(all_cdhits[i])
	sampleid <- unique(datai$sample)
	print(sampleid)
	seqs <- read.delim(paste0(base_dir, "FOR_MINING/", sampleid, "_ALLSAMPLES_WITHSAMPLEID.txt"), header=TRUE)
	seqs <- data.frame(seqs)
	seqs$id <- rownames(seqs)
	seqs$X<- NULL
	datai <- merge(datai, seqs, by.x="gene", by.y="id", all.x=TRUE)
	datai <- merge(datai, reference, by.x="gene", by.y="cdr3_id", all.x="TRUE")
	datai$fullseq <- NULL
	datai$stat <- NULL
	datai$length_cdr3 <- NULL
	## reformat a few of them 
	datai$pathogen[datai$pathogen =="HHV-3 - vaccine"] <- "HHV-3"
	datai$pathogen[datai$pathogen =="HHV-5,HSV  "] <- "HHV-5_HSV"
	datai$pathogen[datai$pathogen =="HHV-2 "] <- "HHV-2"
	datai$pathogen[datai$pathogen =="HHV-3 "] <- "HHV-3"
	datai$pathogen[datai$pathogen =="HHV-5 "] <- "HHV-5"
	datai$pathogen[datai$pathogen =="HHV-4_HHV-2 "] <- "HHV-4_HHV-2"

	#######################
	## Getting the number of references pooled to each cdr3
	seqs_counts <- list.files(paste0(base_dir), full.name=TRUE)
	seqs_counts <- grep(sampleid, seqs_counts, value=TRUE)
	seqs_counts <- grep(".txt", seqs_counts, value=TRUE)
	seqs_counts <- data.frame(rbindlist(sapply(seqs_counts, fread, simplify = FALSE),use.names = TRUE))
	colnames(seqs_counts)[1] <- c("screen_seqs")
	colnames(seqs_counts)[4] <- c("gene")
	
	#######################
	## Creating a look up table of cluster id and all the pathogens detected in that cluster
	## We need to get a dataframe of cluster and pathogen id 
	print("Creating Look up Table")
	datai$pathogen[!datai$gene %like% "cdr3"] <- "HOST"
	s <- datai %>% dplyr::group_by(cluster) %>% mutate(id = paste0(unique(pathogen), collapse = "__")) %>% distinct(cluster, id, .keep_all = TRUE)  
	s <- data.frame(s)
	look_up <- s[, c("cluster", "id")]
	colnames(look_up) <- c("Cluster", "Summary_Pathogen")
	look_up$Sample <- sampleid
	look_up$Summary_Pathogen <- gsub("HOST__", "", look_up$Summary_Pathogen)
	print("Done")

	## Now we want to just get the cdr3s for the screened sequences (excluding reference clusters)
	cdr3s_screen <- datai[!datai$gene %like% "cdr3",]
	cdr3s_screen <- merge(cdr3s_screen, seqs_counts, by=c("gene", "screen_seqs"))		
	## Lets assign the Day to each CDR3.........
	cdr3s_screen$DAY <- NA
	cdr3s_screen$DAY <- str_split_fixed(cdr3s_screen$gene, "\\:", 2)[,1]
	if(unique(cdr3s_screen$sample)%like% "HV"){
		cdr3s_screen$DAY <- str_split_fixed(cdr3s_screen$DAY, "_", 3)[,3]
	} else {
		cdr3s_screen$DAY <- str_split_fixed(cdr3s_screen$DAY, "_", 2)[,2]
	}


	#cdr3s_screen2 <- cdr3s_screen  
	#######################
	## GETTING ISOTYPE COUNTS
	## We want the total reads not unique (take into account if any sequences are clonally expanded). 
	## Take seqid and match it to the total reads for that id 
	files_counts <- list.files(bcr_sequences_iso, full.names=TRUE)
	files_counts <- grep(sampleid, files_counts, value=TRUE)
	files_counts <- grep("_PRODUCTIVE", files_counts, value=TRUE)
	files_counts <- lapply(files_counts, function(x) {read.table(file = x, header = T, sep ="\t")})
	# Combine them
	combined_df <- do.call("rbind", lapply(files_counts, as.data.frame)) 
	combined_df$DAY <- combined_df$Sample
	combined_df$DAY <- gsub("_PRODUCTIVE", "", combined_df$DAY)
	combined_df$DAY <- gsub(paste0(sampleid, "_"), "", combined_df$DAY)
	multis <- cdr3s_screen[cdr3s_screen$cdr3count >1,]
	cdr3s_screen <- merge(cdr3s_screen, combined_df, by=c("seqid", "DAY"), all.x=TRUE)
	
	######################
	## Now for each row we want to get the counts of total isotype rather than cdr3 counts 
	## This can take a bit of time!! :( 
	## multi node lines
	print("Calculating CDR3 Count for aggregated cdr3s")
	for(s in 1:length(rownames(multis))){
		if (s %% 1000 == 0) {print(s)}
		multi1 <- multis[s,]
		gene <- unique(multi1$gene)
		multi_seqs <- data.frame(unlist(str_split(multi1$seqid, "\\,")))
		multi_seqs$DAY <- unique(multi1$DAY)
		colnames(multi_seqs)[1] <- "seqid"
		multi_data <-  merge(combined_df, multi_seqs, by=c("seqid", "DAY"))
		iso_sums <- colSums(multi_data[, 3:14])
		cdr3s_screen[cdr3s_screen$gene==gene,12:23] <- iso_sums
		}
	cdr3s_screen$Sample <- NULL
	print("Done")
	
	######################################
	## duplicate based on ALL ISOTYPE COUNTS!!! Not this wont exactly equal the combined df as some reads removed for being too short!
	## Need to make it numeric so that rep will repeat it 
	cdr3s_screen$ALL <- as.numeric(cdr3s_screen$ALL)
	cdr3s_screen2 <- data.frame(lapply(cdr3s_screen, rep, cdr3s_screen$ALL))
	write.table(cdr3s_screen2, paste0(outputdir, "ISO_USAGE_", sampleid, ".txt"), sep="\t")
	
	## Individual Pathogens (aka sorting for those which are polyreactive!!)
	## We are only doing it for these which we actually found hits for to save time..
	patho_final <- unlist(str_split(all_combinations_pathogen, "_"))
	patho_final <- unlist(str_split(patho_final, ","))
	patho_final <- patho_final[!is.na(patho_final)]
	patho_final <- gsub("HSV  ","HSV", patho_final)
	patho_final <- gsub("HHV-2 ","HHV-2", patho_final)
	
	### We want to speed this up some how 
	### Now we want to subsample and get hits!
	## Lets merge the look up table with the cdr3_screen
	## some pathogen combinations will contain no host sequences in which case we ignore those
	cdr3s_screen2 <- merge(cdr3s_screen, look_up, by.x="cluster", by.y="Cluster")
	cdr3s_screen2 <- as.data.table(cdr3s_screen2)
	
	all_detected <- unique(cdr3s_screen2$Summary_Pathogen)
	all_detected <- unique(unlist(str_split(all_detected, "__")))
	all_detected <- unique(unlist(str_split(all_detected, "_")))
	
	## Check if all are contained within patho final 
	z <- all(all_detected %in% patho_final)
	if(z==TRUE){
		print("All pathogens correctly accounted for")
	 } else {
	  print("WARNING SOMETHING HAS GONE WRONG")
	  break 
	 }
	
	#### Function to do the subsampling
	subsample_function <- function(cdr3s_screen2_day, sub_depth, day_use, x, all_combinations_pathogen){	
			if (x %% 1000 == 0) {print(x)}
			subbed_data <- sample_n(cdr3s_screen2_day, sub_depth)
			clusters <- data.frame(subbed_data[, c("cluster", "Summary_Pathogen")])
			clusters$pathogen <- factor(clusters$Summary_Pathogen, levels=all_combinations_pathogen)
			clusters2 <- clusters %>%group_by(pathogen) %>%summarise(n = n()) %>% tidyr::complete(pathogen, fill = list(n = 0)) %>% mutate(freq = n / sum(n))
			clusters2 <- data.frame(clusters2)
			### Now we want to calculate per individual pathogen - if polyreactive we will count twice once for each pathogen 
			reformat_clusters <- data.frame(matrix(nrow = length(patho_final), ncol = 2))
			colnames(reformat_clusters) <- c("pathogen", "n")
			for(g in 1:length(patho_final)){
				reformat_clusters$pathogen[g] <- patho_final[g]
				reformat_clusters$n[g] <- sum(clusters2$n[clusters2$pathogen %like% patho_final[g]])
			}
			reformat_clusters$freq <- reformat_clusters$n/sub_depth
			reformat_clusters$iteration <- x
			reformat_clusters$DAY <- day_use
			return(reformat_clusters)
	}
	
	print("Starting Subsampling")
	## Actually do the subsampling
	subcounts <- c()
	reps <- 1:10000
	for(f in 1:length(unique(cdr3s_screen2$DAY))){
		day_usex <- unique(cdr3s_screen2$DAY)[f]
		print(paste0("Day ", day_usex))
		cdr3s_screen2_dayx <- cdr3s_screen2[DAY==day_usex,]
		out <- do.call(rbind,mapply(subsample_function, x=reps,MoreArgs = list(cdr3s_screen2_day = cdr3s_screen2_dayx, all_combinations_pathogen=all_combinations_pathogen, day_use=day_usex, sub_depth=sub_depth),SIMPLIFY = FALSE))
		subcounts <- rbind(out, subcounts)
	}	
	subcounts$n <- as.numeric(subcounts$n)
	subcounts$freq <- as.numeric(subcounts$freq)
	mean_hits <- aggregate(list(subcounts$n, subcounts$freq), list(subcounts$pathogen,subcounts$DAY), FUN=mean) 
	colnames(mean_hits) <- c("Pathogen", "Day", "MeanCount", "MeanFreq")
	mean_hits$Sample <- sampleid
	print("DONE")
	
	pdf(paste0(plot_dir, "/SAMPLEHITS_", sampleid, ".pdf"), width=15, height=15)
	plot(ggplot(mean_hits, aes(x=Day, y=log(MeanFreq))) + geom_point() + facet_wrap(~Pathogen, scales="free")+ylab("Log Mean Frequency")+theme_bw())
	dev.off()
	return(mean_hits)
}