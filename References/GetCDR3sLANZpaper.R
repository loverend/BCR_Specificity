## Code to get ebv specific CDR3s from LANZ paper 

setwd('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining')
## cdrs with ebv specificity 
cdr2_info <- read.delim('References/CDR3_EBV_FROM_LANZ_PAPERinfo.txt', sep="\t")
cdr2_info$tidyid <- str_split_fixed(cdr2_info$original_id, "p", 2)[,1]
cdr2_info$barcode_id <- paste0("p",str_split_fixed(cdr2_info$original_id, "p", 2)[,2])


## All the vdj data from the paper
imgt_files <- list.files('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/SpecificityMining/samples_ForRachael', recursive=TRUE, full.name=TRUE)
imgt_files <- grep("AA", imgt_files, value=TRUE)
imgt_files <- grep("CSF", imgt_files, value=TRUE)


## lets read in all the data 
df <- rbindlist(sapply(imgt_files, fread, simplify = FALSE),use.names = TRUE, idcol = "FileName")
df$tidyid <- basename(df$FileName)
df$tidyid <- str_split_fixed(df$tidyid, "_", 2)[,1]

cdr3s <- data.frame(df[, c("tidyid", "barcode_id", "CDR3.IMGT.x")])

#new <- merge(cdr3s, cdr2_info, by=c("tidyid","barcode_id"), all.y=TRUE)
new2 <- merge(cdr3s, cdr2_info, by=c("tidyid","barcode_id"))

## some are missing 
new[is.na(new$CDR3.IMGT.x),] 


new2$fullseq <- new2$CDR3.IMGT.x
new2$CDR3.IMGT.x <- NULL
new2$length_cdr3 <- nchar(new2$fullseq)
min(new2$length_cdr3) ## good we want minimum to be 5!


##########
cdr_final <- new2[, c("fullseq", "length_cdr3", "pathogen", "original_id", "cdr3_id")]
write.table(cdr_final, "EBV_CDR3_final.txt", sep="\t", row.names=FALSE)

cdr3_vec <- cdr_final$fullseq
names(cdr3_vec) <- cdr_final$cdr3_id

write.table(cdr3_vec, "EBV_CDR3_vec.txt", sep="\t", col.names=FALSE)
saveRDS(cdr3_vec, "EBV_CDR3_vec.rds")


# DONE