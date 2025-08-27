#################################################################
# Gene-probe mappings:
#   - target - ensemble
#   - tcga - the gene name is in the format (gene|id), but rows
#       with question marks should be ignored. For example,
#       'A2LD1|87769' has gene name 'A2LD1'
#   - gse37642_1, gse6891_1, gse6891_2 - GPL570
#   - gse37642_2 - GPL96
#   - gse71014 - GPL10558
#   - BEAT_BMA, BEAT_PB - the row names contain the genes
#################################################################

library(tidyverse)
library(rmarkdown)
library(dplyr)
library(tidyr)


target <- readRDS('LEUKEMIA_LIB_data/target.rds')
ensemble <- readRDS('LEUKEMIA_LIB_data/ensemble.rds')

gse6891_1 <- readRDS('LEUKEMIA_LIB_data/GSE6891_1.rds')
load('data_RData/GPL570.RData')

gse6891_2 <- readRDS('LEUKEMIA_LIB_data/GSE6891_2.rds')

gse37642_1 <- readRDS('LEUKEMIA_LIB_data/GSE6891_2.rds')

gse37642_2 <- readRDS('LEUKEMIA_LIB_data/GSE6891_1.rds')
load('data_RData/GPL96.RData')

gse71014 <- readRDS('LEUKEMIA_LIB_data/GSE71014.rds')
load('data_RData/GPL10558.RData')

BEAT_BMA <- readRDS('LEUKEMIA_LIB_data/BEAT_BMA.rds')

BEAT_PB <- readRDS('LEUKEMIA_LIB_data/BEAT_PB.rds')

tcga <- readRDS('LEUKEMIA_LIB_data/TCGA.rds')


# for each dataset, process it to get data for each gene
#   - if there multiple genes in the same dataset --> 
#     use the one with the highest mean expression

# at this point, you have expression data for each gene with no duplicate
# gene names

# generate a frequency table for genes across datasets
# filter the data by removing genes that appear in <80% of the datasets

# Now you have a set of gene expression datasets that all have the
# same genes in them (you may want to make sure the names are consistent --
# in the same order)




encode_probes <- function(dataset, key, symbol_name, ID_name) {

  dataset_matched<-data.frame(dataset$X)
  dataset_matched$gene_ID  <- key[[symbol_name]][match(rownames(dataset_matched), key[[ID_name]])]
  missing_genes <- rownames(dataset_matched[is.na(dataset_matched$gene_ID),])
  rownames(dataset_matched) <- NULL
  dataset_matched <- drop_na(dataset_matched)
  dataset_matched[is.na(dataset_matched$gene_ID),]
  
  dataset_matched %>% 
    mutate(gene_ID=str_split(gene_ID, " */// *")) %>% 
    unnest(gene_ID) ->
    dataset_matched
  
  n_occur <- data.frame(table(dataset_matched$gene_ID))
  nc <- n_occur[n_occur$Freq > 1,]
  max_df <- data.frame()
  dup_names <- as.vector(t(nc$Var1))
  max_df <- dataset_matched %>%
    filter(gene_ID %in% dup_names) %>%
    group_by(gene_ID) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")
  
  dataset_matched <- dataset_matched %>% subset(!gene_ID %in% dup_names)
  return(dataset_matched <- dataset_matched %>% suppressMessages(full_join(max_df)))
}


target_cleaned <- encode_probes(target,ensemble,'SYMBOL','ID')

gse6891_1_cleaned <- encode_probes(gse6891_1,GPL570,'Gene Symbol','ID')
gse6891_2_cleaned <- encode_probes(gse6891_2,GPL570,'Gene Symbol','ID')
gse37642_1_cleaned <- encode_probes(gse37642_1,GPL570,'Gene Symbol','ID')



gse37642_2_cleaned <- encode_probes(gse37642_2,GPL96,'Gene Symbol','ID')



gse71014_cleaned <- encode_probes(gse71014,GPL10558,'Symbol','ID')

tcga_cleaned <- tcga$X[!grepl('?',rownames(tcga$X), fixed = TRUE),]
rownames(tcga_cleaned) <- sub("\\|.*", "", rownames(tcga_cleaned))
tcga_cleaned <- as.data.frame(tcga_cleaned)
tcga_cleaned$gene_ID <- rownames(tcga_cleaned)
rownames(tcga_cleaned) <- NULL

BEAT_BMA$X$gene_ID <- rownames(BEAT_BMA$X)
BEAT_PB$X$gene_ID <- rownames(BEAT_PB$X)
rownames(BEAT_PB$X) <- NULL
rownames(BEAT_BMA$X) <- NULL

BEAT_BMA_cleaned <- BEAT_BMA$X

BEAT_PB_cleaned <- BEAT_PB$X

table(arrange(target_cleaned, gene_ID)$gene_ID)

d1 <- data.frame(gene_ID = arrange(target_cleaned, gene_ID)$gene_ID)
d2 <- data.frame(gene_ID = arrange(gse6891_1_cleaned, gene_ID)$gene_ID)
d3 <- data.frame(gene_ID = arrange(gse6891_2_cleaned, gene_ID)$gene_ID)
d4 <- data.frame(gene_ID = arrange(gse37642_1_cleaned, gene_ID)$gene_ID)
d5 <- data.frame(gene_ID = arrange(gse37642_2_cleaned, gene_ID)$gene_ID)
d6 <- data.frame(gene_ID = arrange(tcga_cleaned, gene_ID)$gene_ID)
d7 <- data.frame(gene_ID = arrange(BEAT_BMA_cleaned, gene_ID)$gene_ID)
d8 <- data.frame(gene_ID = arrange(BEAT_PB_cleaned, gene_ID)$gene_ID)
all_genes <- rbind(d1, d2, d3, d4, d5, d6, d7, d8)

freq_list <- all_genes %>%
  group_by(gene_ID) %>%
  summarise(freq = n()) %>%
  arrange(desc(freq))

#5318 genes left over
gene_overlap <- freq_list[freq_list$freq > 5,]$gene_ID
#2498 genes left over
gene_all <- freq_list[freq_list$freq > 7,]$gene_ID

target_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  target_cleaned

gse6891_1_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  gse6891_1_cleaned

gse6891_2_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  gse6891_2_cleaned

gse37642_1_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  gse37642_1_cleaned

gse37642_2_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  gse37642_2_cleaned

tcga_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  tcga_cleaned

BEAT_BMA_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  BEAT_BMA_cleaned

BEAT_PB_cleaned %>%
  filter(gene_ID %in% gene_overlap) ->
  BEAT_PB_cleaned




target_cleaned <- target_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
gse6891_1_cleaned <- gse6891_1_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
gse6891_2_cleaned <- gse6891_2_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
gse37642_1_cleaned <- gse37642_1_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
gse37642_2_cleaned <- gse37642_2_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
tcga_cleaned <- tcga_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
BEAT_BMA_cleaned <- BEAT_BMA_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
BEAT_PB_cleaned <- BEAT_PB_cleaned %>% column_to_rownames("gene_ID") %>% t() %>% as.data.frame() %>% rownames_to_column("patient")

target_cleaned <- cbind(target_cleaned,risk=target$Y$risk)
gse6891_1_cleaned <- cbind(gse6891_1_cleaned,risk=gse6891_1$Y$risk)
gse6891_2_cleaned <- cbind(gse6891_2_cleaned,risk=gse6891_2$Y$risk)
gse37642_1_cleaned <- cbind(gse37642_1_cleaned,risk=gse37642_1$Y$risk)
gse37642_2_cleaned <- cbind(gse37642_2_cleaned,risk=gse37642_2$Y$risk)
tcga_cleaned <- cbind(tcga_cleaned,risk=tcga$Y$risk)
BEAT_BMA_cleaned <- cbind(BEAT_BMA_cleaned,risk=BEAT_BMA$Y$risk)
BEAT_PB_cleaned <- cbind(BEAT_PB_cleaned,risk=BEAT_PB$Y$risk)

target_cleaned <- droplevels(subset(target_cleaned, risk %in% c("Low", "High")))
gse6891_1_cleaned <- droplevels(subset(gse6891_1_cleaned, risk %in% c("Favorable", "Adverse")))
gse6891_2_cleaned <- droplevels(subset(gse6891_2_cleaned, risk %in% c("Favorable", "Adverse")))
gse37642_1_cleaned <- droplevels(subset(gse37642_1_cleaned, risk %in% c("Favorable", "Adverse")))
gse37642_2_cleaned <- droplevels(subset(gse37642_2_cleaned, risk %in% c("Favorable", "Adverse")))
tcga_cleaned <- droplevels(subset(tcga_cleaned, risk %in% c("favorable", "poor")))
BEAT_BMA_cleaned <- droplevels(subset(BEAT_BMA_cleaned, risk %in% c("Favorable", "Adverse")))
BEAT_PB_cleaned <- droplevels(subset(BEAT_PB_cleaned, risk %in% c("Favorable", "Adverse")))


target_cleaned <- target_cleaned %>% mutate(risk = risk == "High")
gse6891_1_cleaned <- gse6891_1_cleaned %>% mutate(risk = risk == "Adverse")
gse6891_2_cleaned <- gse6891_2_cleaned %>% mutate(risk = risk == "Adverse")
gse37642_1_cleaned <- gse37642_1_cleaned %>% mutate(risk = risk == "Adverse")
gse37642_2_cleaned <- gse37642_2_cleaned %>% mutate(risk = risk == "Adverse")
tcga_cleaned <- tcga_cleaned %>% mutate(risk = risk == "poor")
BEAT_BMA_cleaned <- BEAT_BMA_cleaned %>% mutate(risk = risk == "Adverse")
BEAT_PB_cleaned <- BEAT_PB_cleaned %>% mutate(risk = risk == "Adverse")

combined_all <- bind_rows(
  target_cleaned %>% select(patient, all_of(gene_all), risk),
  gse6891_1_cleaned %>% select(patient, all_of(gene_all), risk),
  gse6891_2_cleaned %>% select(patient, all_of(gene_all), risk),
  gse37642_1_cleaned %>% select(patient, all_of(gene_all), risk),
  gse37642_2_cleaned %>% select(patient, all_of(gene_all), risk),
  tcga_cleaned %>% select(patient, all_of(gene_all), risk),
  BEAT_BMA_cleaned %>% select(patient, all_of(gene_all), risk),
  BEAT_PB_cleaned %>% select(patient, all_of(gene_all), risk)
)

saveRDS(target_cleaned, file = "CleanedDatasets/target_cleaned.rds") 
saveRDS(gse6891_1_cleaned, file = "CleanedDatasets/gse6891_1_cleaned.rds") 
saveRDS(gse6891_2_cleaned, file = "CleanedDatasets/gse6891_2_cleaned.rds") 
saveRDS(gse37642_1_cleaned, file = "CleanedDatasets/gse37642_1_cleaned.rds") 
saveRDS(gse37642_2_cleaned, file = "CleanedDatasets/gse37642_2_cleaned.rds") 
saveRDS(tcga_cleaned, file = "CleanedDatasets/tcga_cleaned.rds") 
saveRDS(BEAT_BMA_cleaned, file = "CleanedDatasets/BEAT_BMA_cleaned.rds") 
saveRDS(BEAT_PB_cleaned, file = "CleanedDatasets/BEAT_PB_cleaned.rds") 
saveRDS(combined_all, file = "CleanedDatasets/combined_all_genes.rds")
