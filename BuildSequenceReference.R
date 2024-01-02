## Formatting the reference of sequences to screen based on my csv file
## Lauren Overend 
## October 2022

#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0


setwd('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining')
## Outputdir from BCR TCR processing pipeline
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
functionality <- "PRODUCTIVE"
runname <- "Screen_Sepsis"

library(seqinr)

#...................................................................
####### Stage1 
####### get pathogen sequences to use for screen  
file = "new_annotated_sequences_nterm_nottrimmed_infectiousorganismsuk.csv"  # This can be edited to decide which sequences to include
p <- as.matrix(read.csv(file, head=T, sep=","))
p <- data.frame(p)
write.table(unique(unique(p$antigen_species)), "Species.txt",  sep="\t", row.names=FALSE)
## We want to remove those we dont want! e.g. annotated with a 1 
p=p[which(p[,"Exclude"]!="1"),]
p <- data.frame(p)
## Lets get just the heavy chain so we can look for cases of duplication 
p <- p[, c("antigen_species", "Hchain_aa")]
## Lets make unique (remove duplicate heavy chains)
d <- unique(p)
f <- table(d$Hchain_aa)
## Lets check no incident of the same species but different length
any(f>1) ## IF TRUE WE KNOW WE HAVE A FEW DUPLICATES where species are different but same sequence (we do!)
n <- data.frame(f[f>1])
dim(n) ## 7 - there are 7 sequences appearing for different organims ## lets exclude these 
## what are the duplicated ones 
e <- d[d$Hchain_aa %in% n$Var1,]
e <- e[order(e$Hchain_aa), ]
seqs_to_remove <- unique(e$Hchain_aa)
## lets remove them 
d <- d[!d$Hchain_aa %in% seqs_to_remove,]
d$length <- nchar(d$Hchain_aa)
d <- d[d$length >5,]
## Lets assign a unique id to each sequence:
d$newID <- paste0("pathogen_", sample(nrow(d)))

seqs_use <- d$Hchain_aa
seqs_usex <- as.list(seqs_use)
names(seqs_use) <- d$newID


## Lets assign a unique id to each sequence



### D is a nice clean dataframe of unique heavy chain sequences. 
#.....................................................................
write.table(d, "Unique_Pathogen_Sequences_meta_info_pretrim.txt",  sep="\t", row.names=FALSE)
saveRDS(seqs_use, "Unique_Pathogen_Sequences_for_screen_pretrim.rds")
write.fasta(seqs_usex, names(seqs_usex), "Unique_Pathogen_Sequences_for_screen_pretrim.fasta", open = "w", nbchar = 60, as.string = FALSE)
#.....................................................................

#....................................................................
## RUN ANARCI!!!!!!!!! to extract CDR3s from protien sequence
#....................................................................

## PART2 ANARCI output 
aligned_file <- read.delim("annotated_cdr3_H.csv", header=TRUE, sep=",")
## IMGT: 
#Complementarity determining regions: CDR1-IMGT: 27 to 38, CDR2-IMGT: 56 to 65 and ***CDR3-IMGT: 105 to 117.******
##Extract CDR3 sequences 
##https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
cdr3_only <- aligned_file %>% select("Id", 'X105':'X117')
rownames(cdr3_only) <- cdr3_only$Id
cdr3_only$Id <- NULL
cdr3_only$fullseq <- apply(cdr3_only, 1, paste, collapse="")
cdr3_only$fullseq <- gsub("-", "", cdr3_only$fullseq)
cdr3_only$length_cdr3 <- nchar(cdr3_only$fullseq)
## Remove cdr3s less than 5aa 
cdr3_only <- cdr3_only[cdr3_only$length_cdr3>5,]
cdr3_only <- cdr3_only[, c("fullseq", "length_cdr3")]
cdr3_only$newID <- rownames(cdr3_only)


## We actually have quite a lot of duplicates in the file e.g. the cdr3 is the same but there is a bit of difference elsewhere (often from trimming)
## Lets read in the info about pathogen so we can concatenate this info 

info <- read.delim("Unique_Pathogen_Sequences_meta_info_pretrim.txt", sep="\t")
info2 <- merge(cdr3_only, info, by="newID")
info2 <- info2[, c("newID", "fullseq", "antigen_species", "length_cdr3")]
info3 <- info2[, c("fullseq", "length_cdr3")]
## Now we have only unique CDR3s 
unique_cdr3 <- unique(info3)
unique_cdr3$pathogen <- NA
unique_cdr3$originalID <- NA

for(i in 1:length(unique_cdr3[,1])){
	cdr3_use <- unique_cdr3$fullseq[i]
	datax <- info2[info2$fullseq==cdr3_use,]
	pathogens <- unique(datax$antigen_species)
	pathogens <- paste(pathogens,collapse=",")
	seqids <- unique(datax$newID)
	seqids <- paste(seqids,collapse=",")
	unique_cdr3$pathogen[i] <- pathogens
	unique_cdr3$originalID[i] <- seqids
}

unique_cdr3$cdr3_id <- paste0("cdr3_", sample(nrow(unique_cdr3)))

### get the unique cdr3s 
cdr3 <- unique_cdr3$fullseq
names(cdr3) <- unique_cdr3$cdr3_id


saveRDS(cdr3, "Unique_Pathogen_Sequences_for_screen_POST_trim.rds")
write.table(unique_cdr3, "Unique_CDR3_info.txt",  sep="\t", row.names=FALSE)

#....................................................................
## DONE - reference tidied, trimmed, made unique and saved. 
#....................................................................

## Lets add on the references from the LANZ paper 

cdr3 <- readRDS('References/Unique_Pathogen_Sequences_for_screen_POST_trim.rds')
cdr3_ebv <- readRDS("References/EBV_CDR3_vec.rds")
all_ref <- c(cdr3, cdr3_ebv)

saveRDS(all_ref, "HUMAN_PATHOGEN_UK_REFERENCE.rds")

info1 <- read.delim("References/Unique_CDR3_info.txt")
info2 <- read.delim("References/EBV_CDR3_final.txt")
colnames(info2)[4] <- "originalID"
final <- rbind(info1, info2)
write.table(final, "HUMAN_PATHOGEN_UK_REFERENCE_INFO.txt", sep="\t", row.names=FALSE, quote=FALSE)

### DONE

