## Analysis of specificity mining hits 
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk


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

#............................................................................................................
### ANALYSIS FOR CDHIT 

outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/CDHIT_PERINDIVIDUAL/'
results <- list.files(outputdir, full.names=TRUE)
results <- grep("screen.clstr", results, value=TRUE)
reference <- read.delim('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining/HUMAN_PATHOGEN_UK_REFERENCE_INFO.txt')
plot_dir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/CDHIT_PERINDIVIDUAL/PLOTS'
base_dir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/'
bcr_sequences_iso <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/SummarySequences'
metadata <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt', header=TRUE)
metadata <- metadata[, c( "SampleID_alternative", "alternative_barcode", "Mortality.Classification", "Sex", "Age")]
metadata$Mortality2 <- metadata$Mortality.Classification
metadata$Mortality2[metadata$Mortality2=="MID.DEATH"] <- "LATE.DEATH"
setwd('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining/')
## Lets remove the technical 
results <- results[!results %like% "JR1795"]

#...................................................................................................	
## Part 1:
## Get the CDHIT results and reformat into a dataframe - saving for each individual.
source("FormatCDHITResults.R")
chhits <- c()
for(z in 1:length(results)){
	cdhits <- format_cdhit(results, z)
	chhits <- rbind(cdhits, chhits)
}
## Simplify some of the pathogen ids for later!!!
chhits$pathogen[chhits$pathogen =="HHV-3 - vaccine"] <- "HHV-3"
chhits$pathogen[chhits$pathogen =="HHV-5,HSV  "] <- "HHV-5_HSV"
chhits$pathogen[chhits$pathogen =="HHV-3 - vaccine"] <- "HHV-3"
chhits$pathogen[chhits$pathogen =="HHV-5,HSV  "] <- "HHV-5_HSV"
chhits$pathogen[chhits$pathogen =="HHV-2 "] <- "HHV-2"
chhits$pathogen[chhits$pathogen =="HHV-3 "] <- "HHV-3"
chhits$pathogen[chhits$pathogen =="HHV-5 "] <- "HHV-5"
chhits$pathogen[chhits$pathogen =="HHV-4_HHV-2 "] <- "HHV-4_HHV-2"
write.table(chhits, paste0(outputdir, "/SUMMARY_SEPSIS_CDHits.txt"), sep="\t")
print("Done")
###### Load so we dont need to keep rerunning 
chhits <- read.delim( paste0(outputdir, "/SUMMARY_SEPSIS_CDHits.txt"), sep="\t")
all_combinations_pathogen <- c(unique(chhits$pathogen), "HOST")
all_combinations_pathogen <- all_combinations_pathogen[!is.na(all_combinations_pathogen)]

#...................................................................................................	
## Now lets get the total depth which we mined for these samples (per day!!)
## Extract the minimum read depth for downstream subsampling. 
outputdir_readepth <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/CDR3_AA/FOR_MINING'
all_rds1 <- list.files(outputdir_readepth, full.name=TRUE)
all_rds1 <- grep("_WITHSAMPLEID.txt", all_rds1, value=TRUE)
all_rds1 <- grep("PRODUCTIVE", all_rds1, value=TRUE)
read_depth <- c()
for(i in 1:length(all_rds1)){
	sampleid <- basename(all_rds1)[i]
	sampleid <- gsub("_ALLSAMPLES_WITHSAMPLEID.txt", "", sampleid)
	files_all <- read.delim(all_rds1[i])
	readdepth <- length(files_all$screen_seqs[files_all$X=="screen"])
	no_samples <- rownames(files_all)[files_all$X=="screen"]
	no_samples <- length(unique(str_split_fixed( no_samples,"\\:", 2)[,1]))
	readdepth <- c(sampleid, readdepth, no_samples)
	read_depth <- rbind(read_depth, readdepth)
}
read_depth <- data.frame(read_depth)
colnames(read_depth) <- c("Sample", "TotalSeqs", "NoSamples")

#...................................................................................................	
## Now we want to subsample results 
read_depth$TotalSeqs <- as.numeric(read_depth$TotalSeqs)
min_depth <- min(read_depth$TotalSeqs)
sub_depth <- 0.8*as.numeric(min_depth)
sub_depth <- round(sub_depth)
read_depth$DISEASE <- "SEPSIS"
read_depth$DISEASE[read_depth$Sample %like% "HV"] <- "HEALTH"
## Histogram of read depths 
pdf(paste0(plot_dir, "/ReadDepth.pdf"), width=6, height=4)
p<-ggplot(read_depth, aes(x=TotalSeqs, fill=DISEASE, colour=DISEASE)) + geom_histogram(position="identity", alpha=0.5) +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")
plot(p)
dev.off()
## Healthy Controls have much higher read depth so higher odds of getting a hit! ## We need to subsample the matrix 

#...................................................................................................
## PART 2, SUBSAMPLE TO FIXED DEPTH AND CALCULATE HITS 
all_cdhits <- list.files(outputdir, full.name=TRUE)
all_cdhits <- grep("Hits", all_cdhits, value=TRUE)
## Custom code to subsample! 
## This is going to take a very long time!!! :(
source("SubsampleCDHITS.R")
pathogen_hits <- c() 
for(i in 1:length(all_cdhits)){
	hits <- subsample_cdhits(all_cdhits, i, base_dir, bcr_sequences_iso, sub_depth, plot_dir)
	pathogen_hits <- rbind(hits, pathogen_hits)
}
write.table(pathogen_hits, paste0(base_dir, "/SEPSIS_CDHits_SUBSAMPLED_10k.txt"), sep="\t")
print("Done")
pathogen_hits <- read.delim("SEPSIS_CDHits_SUBSAMPLED_10k.txt")
pathogen_hits$SampleID <- paste0(pathogen_hits$Sample, "_", pathogen_hits$Day)

##...............................................


pathogen_hits2 <- merge(pathogen_hits, metadata, by.x="SampleID", by.y="SampleID_alternative", all.x=TRUE)
pathogen_hits2$Mortality2[pathogen_hits2$Sample %like% "HV"] <- "HEALTH"
pathogen_hits2$DISEASE <- "SEPSIS"
pathogen_hits2$DISEASE[pathogen_hits2$Sample %like% "HV"] <- "HEALTH"


pdf(paste0(plot_dir, "/PathogenHits.pdf"), width=15, height=15)
p<-ggplot(pathogen_hits2, aes(x=as.factor(Day), y=log(MeanFreq), fill=DISEASE)) + geom_boxplot() +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)

p<-ggplot(pathogen_hits2, aes(x=Day, y=log(MeanFreq), colour=DISEASE, fill=DISEASE)) + geom_point() +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_smooth(method='lm')
plot(p)

p<-ggplot(pathogen_hits2, aes(x=as.factor(Day), y=MeanFreq, fill=DISEASE)) + geom_boxplot() +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)

p<-ggplot(pathogen_hits2, aes(x=Day, y=MeanFreq, colour=DISEASE, fill=DISEASE)) + geom_point() +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_smooth(method='lm')
plot(p)

dev.off()



p<-ggplot(pathogen_hits2, aes(x=as.factor(Day), y=log(MeanFreq), fill=DISEASE)) + geom_violin(trim=FALSE) +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free") + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)
dev.off()


pdf(paste0(plot_dir, "/PathogenHits_Count.pdf"), width=15, height=15)
p<-ggplot(pathogen_hits2, aes(x=as.factor(Day), y=(MeanCount), fill=DISEASE)) + geom_boxplot() +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)
p<-ggplot(pathogen_hits2, aes(x=as.factor(Day), y=(MeanCount), fill=DISEASE)) + geom_violin(trim=FALSE) +theme_bw() + theme(legend.position="bottom")+xlab("Total Sequences per Individual")+facet_wrap(~Pathogen, scales="free") + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)
dev.off()










#.................................	
## now lets merge depths with hits 

read_depth <- read_depth[!read_depth$Sample %like% "JR1795",]
final <- merge(chhits, read_depth, by.x="sample", by.y="Sample", all=TRUE)

summary_df <- c()
## Lets calculate hits per total read depth and make better format 
for(i in 1:length(unique(final$sample))){
	subs <- final[final$sample==unique(final$sample)[i],]
	if(dim(subs)[1]>1){
		for(x in 1:length(unique(subs$cluster))){
			clus <- subs[subs$cluster==unique(subs$cluster)[x],]
			clusid <- paste0(unique(clus$sample), "_", unique(clus$cluster))
			no_hits <- length(clus$gene[!clus$gene %like% "cdr3"])
			pathogen <- unique(clus$pathogen)
			pathogen <- pathogen[!is.na(pathogen)] 
			## if the cdr3 is polyreactive we want a count for both pathogens
			if(length(pathogen)>1){
				#print(i)
				#print(x)
				total_reads <- unique(clus$TotalSeqs)
				total_samples <- unique(clus$NoSamples)
				sampleid <- unique(clus$sample)
				for(y in 1:length(pathogen)){
					rowx <- c(sampleid, paste0(clusid, ":", y), pathogen[y], no_hits, total_reads, total_samples)
					#print(rowx)
					summary_df <- rbind(summary_df, rowx)
				}
			}
			### If its not polyreactive 
			if(length(pathogen)==1){
				pathogen<- paste(pathogen, collapse=" + ")
				total_reads <- unique(clus$TotalSeqs)
				total_samples <- unique(clus$NoSamples)
				sampleid <- unique(clus$sample)
				rowx <- c(sampleid, clusid, pathogen, no_hits, total_reads, total_samples)
				summary_df <- rbind(summary_df, rowx)
			}
			}
		} else {
		rowx <- c(subs$sample, NA, "NoHits", 1, subs$TotalSeqs, subs$NoSamples)
		summary_df <- rbind(summary_df, rowx)
		}
	}
summary_df <- data.frame(summary_df)
colnames(summary_df) <- c("Individual", "ClusterID", "Pathogen", "Count", "TotalReads", "NoSamples")
summary_df$Count <- as.numeric(summary_df$Count)
summary_df$Individual <- as.factor(summary_df$Individual)
summary_df$Pathogen <- as.factor(summary_df$Pathogen)

## If there was no hit this is a count of 0 hence we want to create a row for this in the df
summary_df2 <- summary_df %>% group_by(Individual, Pathogen, .drop = FALSE) %>% summarise(count=sum(Count))
summary_df2 <- data.frame(summary_df2)
summary_df3 <- summary_df2 
## Remove level of no hits
summary_df2 <- summary_df2[!summary_df2$Pathogen=="NoHits",]
summary_df2$Pathogen <- as.factor(as.character(summary_df2$Pathogen))
summary_df2 <- merge(summary_df2, read_depth, by.x="Individual", by.y="Sample")
summary_df2$TotalSeqs <- as.numeric(summary_df2$TotalSeqs)
summary_df2$Per <- (summary_df2$count/summary_df2$TotalSeqs*100)

metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
metadata <- read.delim(metadata, sep="\t", header=TRUE)
metadata <- metadata[, c( "SampleID_alternative", "alternative_barcode","Mortality.Classification","Sex", "Age")]
metadata$Mortality2 <- metadata$Mortality.Classification
metadata$Mortality2[metadata$Mortality2=="MID.DEATH"] <- "LATE.DEATH"

summary_df2 <- merge(summary_df2, metadata, by.x="Individual", by.y="alternative_barcode", all.x=TRUE)
summary_df2$DISEASE <- "SEPSIS"
summary_df2$DISEASE[summary_df2$Individual %like% "HV"] <- "HEALTH"
summary_df2$Mortality2[summary_df2$DISEASE %like% "HEALTH"] <- "HEALTH"

pdf(paste0(plot_dir, "/CDHIT_MATCHES_LEO_SEPSIS.pdf"), width=15, height=15)
ggplot(summary_df2, aes(x=Mortality2, y=Per, fill=DISEASE))+geom_boxplot() + facet_wrap(~Pathogen, scales="free_y")+theme_bw()+ylab("CDR3 Hits %")+ggtitle("All Samples")+ stat_compare_means(aes(label = ..p.signif..), method="anova")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(summary_df2[summary_df2$Per !=0,], aes(x=Mortality2, y=Per, fill=DISEASE))+geom_boxplot() + facet_wrap(~Pathogen, scales="free_y")+theme_bw()+ylab("CDR3 Hits %")+ggtitle("'Positive' Samples")+ stat_compare_means(aes(label = ..p.signif..), method="anova")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(summary_df2, aes(x=DISEASE, y=Per, fill=DISEASE))+geom_boxplot() + facet_wrap(~Pathogen, scales="free_y")+theme_bw()+ylab("CDR3 Hits %")+ggtitle("All Samples")+ stat_compare_means(aes(label = ..p.signif..), method="wilcox")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(summary_df2[summary_df2$Per !=0,], aes(x=DISEASE, y=Per, fill=DISEASE))+geom_boxplot() + facet_wrap(~Pathogen, scales="free_y")+theme_bw()+ylab("CDR3 Hits %")+ggtitle("'Positive' Samples")+ stat_compare_means(aes(label = ..p.signif..), method="wilcox")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##........................................................
## Calculate proportion of individuals which are positive in Health vs Disease 
summary_df3 <- merge(summary_df3, metadata, by.x="Individual", by.y="alternative_barcode", all.x=TRUE)
summary_df3$DISEASE <- "SEPSIS"
summary_df3$DISEASE[summary_df3$Individual %like% "HV"] <- "HEALTH"
summary_df3$Mortality2[summary_df3$DISEASE %like% "HEALTH"] <- "HEALTH"
summary_df3$hit <- 0
summary_df3$hit[summary_df3$count>0] <- 1

proper <- summary_df3 %>% group_by(DISEASE, Pathogen, hit) %>%summarise(n = n()) %>%mutate(freq = n / sum(n))
proper <- data.frame(proper)
props_positivex1 <- proper[proper$hit==1,]

proper <- summary_df3 %>% group_by(Mortality2, Pathogen, hit) %>%summarise(n = n()) %>%mutate(freq = n / sum(n))
proper <- data.frame(proper)
props_positivex2 <- proper[proper$hit==1,]

pdf("CDHIT_MATCHES_LEO_SEPSIS_PROP_POSITIVE.pdf", width=10, height=10)
ggplot(props_positivex1, aes(x=DISEASE, y=freq, fill=DISEASE))+geom_bar(stat="identity") + facet_wrap(~Pathogen)+theme_bw()+ylab("Proportion of Samples with at least one CDR3 Hit")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(props_positivex2, aes(x=Mortality2, y=freq, fill=Mortality2))+geom_bar(stat="identity") + facet_wrap(~Pathogen)+theme_bw()+ylab("Proportion of Samples with at least one CDR3 Hit")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()