

#Load libraries
library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("kableExtra") 


fastq_dir <- "//vortexfs1/home/bfowler/mvco_fastq/"  # fastq directory
filtered_dir <- "//vortexfs1/home/bfowler/Pipeline2022/fastq_filtered/"  # fastq filtered
qual_dir <- "//vortexfs1/home/bfowler/Pipeline2022/qual_pdf/"  # qual pdf
dada2_dir <- "//vortexfs1/home/bfowler/Pipeline2022/dada2/"  # dada2 results
blast_dir <- "//vortexfs1/home/bfowler/Pipeline2022/blast/"  # blast2 results
database_dir <- "//vortexfs1/home/bfowler/databases/"  # databases


PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", 
                    "Genus", "Species")

# get a list of all fastq files in the ngs directory and separate R1 and R2
fns <- sort(list.files(fastq_dir, full.names = TRUE, recursive = TRUE))
fns_batch1 <- fns[str_detect(basename(fns), "SRR")]

#this is actually two batches separated in teh worst fuckign way 
batch1list <- c("SRR8175836",
"SRR8175837",
"SRR8175838",
"SRR8175839",
"SRR8175840",
"SRR8175841",
"SRR8175842",
"SRR8175843",
"SRR8175844",
"SRR8175845",
"SRR8175846",
"SRR8175847",
"SRR8175848",
"SRR8175849",
"SRR8175850",
"SRR8175852",
"SRR8175853",
"SRR8175872",
"SRR8175873",
"SRR8175878",
"SRR8175879",
"SRR8175880",
"SRR8175881",
"SRR8175882",
"SRR8175883",
"SRR8175884",
"SRR8175885")


fns_batch2 <-  fns_batch1[str_detect(basename(fns_batch1), paste(batch1list, collapse = "|"), negate = TRUE)]
fns_batch1 <- fns[str_detect(basename(fns), paste(batch1list, collapse = "|"))]


fns_batch3 <- fns[str_detect(basename(fns), "TH")]
fns_batch4 <- fns[str_detect(basename(fns), "-MV")]


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq==
sample1.names <- str_split(basename(fns_batch1), pattern = ".fastq", simplify = TRUE)
sample1.names <- sample1.names[, 1]

sample2.names <- str_split(basename(fns_batch2), pattern = ".fastq", simplify = TRUE)
sample2.names <- sample2.names[, 1]

sample3.names <- str_split(basename(fns_batch3), pattern = ".fastq", simplify = TRUE)
sample3.names <- sample3.names[, 1]

sample4.names <- str_split(basename(fns_batch4), pattern = ".fastq", simplify = TRUE)
sample4.names <- sample4.names[, 1]

# create an empty data frame
df_1 <- data.frame()
df_2 <- data.frame()
df_3 <- data.frame()
df_4 <- data.frame()


# loop through all the R1 files (no need to go through R2 which should be
# the same)

for (i in 1:length(fns_batch1)) {
  
  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_batch1[i])
  
  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_batch1[i]))  #modified
  
  # add one line to data frame
  df_1 <- bind_rows(df_1, df_one_row)
}

for (i in 1:length(fns_batch2)) {
  
  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_batch2[i])
  
  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_batch2[i]))  #modified
  
  # add one line to data frame
  df_2 <- bind_rows(df_2, df_one_row)
}

for (i in 1:length(fns_batch3)) {
  
  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_batch3[i])
  
  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_batch3[i]))  #modified
  
  # add one line to data frame
  df_3 <- bind_rows(df_3, df_one_row)
}


for (i in 1:length(fns_batch4)) {
  
  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_batch4[i])
  
  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_batch4[i]))  #modified
  
  # add one line to data frame
  df_4 <- bind_rows(df_4, df_one_row)
}


write.table(df_1, file = 'n_seq_batch1.txt', sep='\t', row.names = FALSE, na='', quote=FALSE)
write.table(df_2, file = 'n_seq_batch2.txt', sep='\t', row.names = FALSE, na='', quote=FALSE)
write.table(df_3, file = 'n_seq_batch3.txt', sep='\t', row.names = FALSE, na='', quote=FALSE)
write.table(df_4, file = 'n_seq_batch4.txt', sep='\t', row.names = FALSE, na='', quote=FALSE)


# plot the histogram with number of sequences
#pdf("plot1.pdf")
#myplot <- ggplot(df, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 50000) + xlim(50000, 1000000)
#print(myplot)
#qdev.off()


#for (i in 1:length(fns)) {
# # Use dada2 function to plot quality
#  p1 <- plotQualityProfile(fns[i])
  
#  # Only plot on screen for first 2 files
#  if (i <= 2) {
#    print(p1)
#  }
#  
  # save the file as a pdf file (uncomment to execute)
#  p1_file <- paste0(qual_dir, basename(fns[i]), ".qual.pdf")
 # ggsave(plot = p1, filename = p1_file, device = "pdf", width = 15, height = 15, 
 #        scale = 1, units = "cm")
#}



filt_R1_batch1 <- str_c(filtered_dir, sample1.names, "_R1_filt.fastq")
filt_R1_batch2 <- str_c(filtered_dir, sample2.names, "_R1_filt.fastq")
filt_R1_batch3 <- str_c(filtered_dir, sample3.names, "_R1_filt.fastq")
filt_R1_batch4 <- str_c(filtered_dir, sample4.names, "_R1_filt.fastq")


#trim without Reverse reads
out_1 <- filterAndTrim(fns_batch1, filt_R1_batch1, truncLen = 190,  trimLeft = 17, maxN = 0, maxEE = 2, truncQ = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
out_2 <- filterAndTrim(fns_batch2, filt_R1_batch2, truncLen = 190,  trimLeft = 17, maxN = 0, maxEE = 2, truncQ = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
out_3 <- filterAndTrim(fns_batch3, filt_R1_batch3, truncLen = 190,  trimLeft = 17, maxN = 0, maxEE = 2, truncQ = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
out_4 <- filterAndTrim(fns_batch4, filt_R1_batch4, truncLen = 210,  trimLeft = 17, maxN = 0, maxEE = 2, truncQ = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)


err_1 <- learnErrors(filt_R1_batch1, multithread = FALSE)
err_2 <- learnErrors(filt_R1_batch2, multithread = FALSE)
err_3 <- learnErrors(filt_R1_batch3, multithread = FALSE)
err_4 <- learnErrors(filt_R1_batch4, multithread = FALSE)


#For now just doing forward reads

derep_R1_batch1 <- derepFastq(filt_R1_batch1, verbose = FALSE)
names(derep_R1_batch1) <- sample1.names
derep_R1_batch2 <- derepFastq(filt_R1_batch2, verbose = FALSE)
names(derep_R1_batch2) <- sample2.names
derep_R1_batch3 <- derepFastq(filt_R1_batch3, verbose = FALSE)
names(derep_R1_batch3) <- sample3.names
derep_R1_batch4 <- derepFastq(filt_R1_batch4, verbose = FALSE)
names(derep_R1_batch4) <- sample4.names


dada_R1_batch1 <- dada(derep_R1_batch1, err = err_1, multithread = FALSE, pool = FALSE)
dada_R1_batch2 <- dada(derep_R1_batch2, err = err_2, multithread = FALSE, pool = FALSE)
dada_R1_batch3 <- dada(derep_R1_batch3, err = err_3, multithread = FALSE, pool = FALSE)
dada_R1_batch4 <- dada(derep_R1_batch4, err = err_4, multithread = FALSE, pool = FALSE)


seqtab_1 <- makeSequenceTable(dada_R1_batch1) #using derep instead of merge 
seqtab_2 <- makeSequenceTable(dada_R1_batch2) #using derep instead of merge 
seqtab_3 <- makeSequenceTable(dada_R1_batch3) #using derep instead of merge 
seqtab_4 <- makeSequenceTable(dada_R1_batch4) #using derep instead of merge 

# Make a transposed of the seqtab to make it be similar to mothur database
t_seqtab_1 <- t(seqtab_1)
t_seqtab_2 <- t(seqtab_2)
t_seqtab_3 <- t(seqtab_3)
t_seqtab_4 <- t(seqtab_4)

# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))

seqtab1.nochim <- removeBimeraDenovo(seqtab_1, method = "consensus", multithread = FALSE, 
                                    verbose = TRUE)
seqtab2.nochim <- removeBimeraDenovo(seqtab_2, method = "consensus", multithread = FALSE, 
                                     verbose = TRUE)
seqtab3.nochim <- removeBimeraDenovo(seqtab_3, method = "consensus", multithread = FALSE, 
                                     verbose = TRUE)
seqtab4.nochim <- removeBimeraDenovo(seqtab_4, method = "consensus", multithread = FALSE, 
                                     verbose = TRUE)

# Compute % of non chimeras
paste0("% of non chimeras : ", sum(seqtab1.nochim)/sum(seqtab_1) * 100)
paste0("% of non chimeras : ", sum(seqtab2.nochim)/sum(seqtab_2) * 100)
paste0("% of non chimeras : ", sum(seqtab3.nochim)/sum(seqtab_3) * 100)
paste0("% of non chimeras : ", sum(seqtab4.nochim)/sum(seqtab_4) * 100)


# define a function. changed so no merging 
getN <- function(x) sum(getUniques(x))

track_1 <- cbind(out_1, sapply(dada_R1_batch1, getN), rowSums(seqtab_1), 
               rowSums(seqtab1.nochim))
track_2 <- cbind(out_2, sapply(dada_R1_batch2, getN), rowSums(seqtab_2), 
                 rowSums(seqtab2.nochim))
track_3 <- cbind(out_3, sapply(dada_R1_batch3, getN), rowSums(seqtab_3), 
                 rowSums(seqtab3.nochim))
track_4 <- cbind(out_4, sapply(dada_R1_batch4, getN), rowSums(seqtab_4), 
                 rowSums(seqtab4.nochim))

colnames(track_1) <- c("input", "filtered", "denoised", "tabled", "nonchim")
colnames(track_2) <- c("input", "filtered", "denoised", "tabled", "nonchim")
colnames(track_3) <- c("input", "filtered", "denoised", "tabled", "nonchim")
colnames(track_4) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track_1) <- sample1.names
rownames(track_2) <- sample2.names
rownames(track_3) <- sample3.names
rownames(track_4) <- sample4.names


write_tsv(data.frame(track_1), str_c(dada2_dir, "read_numbers_dada2_batch1.tsv"))
write_tsv(data.frame(track_2), str_c(dada2_dir, "read_numbers_dada2_batch2.tsv"))
write_tsv(data.frame(track_3), str_c(dada2_dir, "read_numbers_dada2_batch3.tsv"))
write_tsv(data.frame(track_4), str_c(dada2_dir, "read_numbers_dada2_batch4.tsv"))


seqtab1.nochim_trans <- as.data.frame(t(seqtab1.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d", 
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
seqtab2.nochim_trans <- as.data.frame(t(seqtab2.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d", 
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
seqtab3.nochim_trans <- as.data.frame(t(seqtab3.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d", 
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
seqtab4.nochim_trans <- as.data.frame(t(seqtab4.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d", 
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

df_1 <- seqtab1.nochim_trans
df_2 <- seqtab2.nochim_trans
df_3 <- seqtab3.nochim_trans
df_4 <- seqtab4.nochim_trans

seq_out_batch1 <- Biostrings::DNAStringSet(df_1$sequence)
seq_out_batch2 <- Biostrings::DNAStringSet(df_2$sequence)
seq_out_batch3 <- Biostrings::DNAStringSet(df_3$sequence)
seq_out_batch4 <- Biostrings::DNAStringSet(df_4$sequence)

names(seq_out_batch1) <- df_1$OTUNumber
names(seq_out_batch2) <- df_2$OTUNumber
names(seq_out_batch3) <- df_3$OTUNumber
names(seq_out_batch4) <- df_4$OTUNumber




Biostrings::writeXStringSet(seq_out_batch1, str_c(dada2_dir, "MVCO_ASV_no_taxo_batch1.fasta"), 
                            compress = FALSE, width = 20000)
Biostrings::writeXStringSet(seq_out_batch2, str_c(dada2_dir, "MVCO_ASV_no_taxo_batch2.fasta"), 
                            compress = FALSE, width = 20000)
Biostrings::writeXStringSet(seq_out_batch3, str_c(dada2_dir, "MVCO_ASV_no_taxo_batch3.fasta"), 
                            compress = FALSE, width = 20000)
Biostrings::writeXStringSet(seq_out_batch4, str_c(dada2_dir, "MVCO_ASV_no_taxo_batch4.fasta"), 
                            compress = FALSE, width = 20000)



pr2_file <- paste0(database_dir, "pr2_version_4.14.0_SSU_dada2.fasta.gz")

save.image(file = "Outputs4batches_all.Rdata")


mergtab.nochim <- mergeSequenceTables(seqtab1.nochim, seqtab2.nochim, seqtab3.nochim, seqtab4.nochim)
mergtab.nochim_trans <- as.data.frame(t(mergtab.nochim)) %>% rownames_to_column(var = "sequence") %>% 
    rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d", 
                                                                    OTUNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

taxa_4batch <- assignTaxonomy(mergtab.nochim, refFasta = pr2_file, taxLevels = PR2_tax_levels, 
                  minBoot = 0, outputBootstraps = TRUE, verbose = TRUE)

saveRDS(taxa_4batch, str_c(dada2_dir, "MVCO_4batch.taxa.rds"))


write_tsv(as.tibble(taxa_4batch$tax), path = str_c(dada2_dir, "taxa4batch.txt"))
write_tsv(as.tibble(taxa_4batch$boot), path = str_c(dada2_dir, "taxa4_boot.txt"))
write_tsv(as.tibble(mergtab.nochim), path = str_c(dada2_dir, "taxa4tab.txt"))

taxa4_tax <- as.data.frame(taxa_4batch$tax)
taxa4_boot <- as.data.frame(taxa_4batch$boot) %>% rename_all(funs(str_c(., "_boot")))
mergtab.nochim_trans <- taxa4_tax %>% bind_cols(taxa4_boot) %>% bind_cols(mergtab.nochim_trans)


write_tsv(as.tibble(mergtab.nochim_trans), file = str_c(dada2_dir, "mergtab4_trans.txt"))


save.image(file = "Outputs4batches_all.Rdata")
