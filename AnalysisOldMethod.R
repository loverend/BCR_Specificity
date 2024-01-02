outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/RESULTS/'
results <- list.files(outputdir, full.names=TRUE)
reference <- read.delim('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining/Unique_CDR3_info.txt')
plot_dir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/Plots'

all_data <- c()
for(i in 1:length(results)){

	## Check for results
	info = file.info(results[i])
	if(info$size==0){
		next 
	}
	
	data_use <- read.delim(results[i])

	data_usex <- merge(data_use, reference, by.x="X.ref_seq", by.y="fullseq")
	sampleid <- basename(results[i])
	sampleid <- gsub("_PRODUCTIVE_CDR3AA_matched_cdr3s.txt", "", sampleid)
	
	## Pathogens with 0 mismatch 
	data_use0 <- data_usex[data_usex$difference ==0,]
	hits0 <- data_use0 %>% group_by(pathogen) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
	hits0 <- data.frame(hits0)
	if(dim(hits0)[1]==0){
		hits0<- hits0 %>% add_column("mismatch" = NA) 
	} else {
		hits0$mismatch <- 0
	}

	## Pathogens with less than 1 mismatch 
	data_use1 <- data_usex[data_usex$difference ==1,]
	hits1 <- data_use0 %>% group_by(pathogen) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
	hits1 <- data.frame(hits1)
	if(dim(hits1)[1]==0){
		hits1<- hits1 %>% add_column("mismatch" = NA) 
	} else {
		hits1$mismatch <- 1
	}
	
	## Pathogens with less than 2 mismatch 
	data_use2 <- data_usex[data_usex$difference ==2,]
	hits2 <- data_use2 %>% group_by(pathogen) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
	hits2 <- data.frame(hits2)
	if(dim(hits2)[1]==0){
		hits2<- hits2 %>% add_column("mismatch" = NA) 
	} else {
		hits2$mismatch <- 2
	}
	
	## Pathogens with less than 1 mismatch 
	data_use3 <- data_usex[data_usex$difference ==3,]
	hits3 <- data_use3 %>% group_by(pathogen) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
	hits3 <- data.frame(hits3)
	if(dim(hits3)[1]==0){
		hits3<- hits3 %>% add_column("mismatch" = NA) 
	} else {
		hits3$mismatch <- 3
	}
	
	all_hits <- rbind(hits0, hits1, hits2, hits3)
	
	## lets do a plot for this sample 
	if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
	
	all_hits$percent <- round(all_hits$freq*100)
	all_hits$sample <- sampleid
	all_hits$mismatch <- as.factor(all_hits$mismatch)
	pdf(paste0(plot_dir, "/HITS_", sampleid, ".pdf"), height=7, width=9)
	x <- ggplot(all_hits, aes(fill=pathogen, y=freq, x=mismatch)) + geom_bar(position="stack", stat="identity", colour="black") +guides(fill=guide_legend(ncol=2))+theme_classic()+xlab("CDR3 mismtach") + ylab("Proportion of hits")  +ggtitle(sampleid) 
	plot(x)
	dev.off()
	
	all_data <- rbind(all_hits, all_data)
}
	