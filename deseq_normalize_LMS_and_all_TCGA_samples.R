# analysis of all TCGA samples and only leio with survival curves
### we aim to correlate or 33 genes signature with survival of patients with these genes de-regulated
{
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
}

# Case ID is the link between the rnaseq file and the data in clinical 
samplesheet_TCGA_LMS <- 'manual_downloaded_samples/gdc_sample_sheet.2025-09-19.tsv'
samplesheet_TCGA_ALL <- 'all_TCGA/gdc_sample_sheet.2025-10-13.tsv'

files_TCGA_LMS <- list.files('manual_downloaded_samples/rnaseq', recursive = T, pattern = '.tsv', full.names = T)
files_TCGA_ALL <- list.files('all_TCGA/data', recursive = T, pattern = '.augmented_star_gene_counts.tsv$', full.names = T)


##########
# TCGA LMS
##########

samplesheet_TCGA_LMS <- read_tsv(samplesheet_TCGA_LMS) %>% 
  dplyr::select(`File Name`, `Case ID`)

tab <- tibble()

for (data in files) {
  df <- read_delim(data, comment = '#', skip = 6, col_names = F)
  df <- df |> 
    dplyr::select(X2, X4) |> 
    dplyr::rename('symbol' = 'X2', 'Raw_counts' = 'X4') |> 
    mutate(sample = basename(data)) |> 
    dplyr::filter(!duplicated(symbol))
  
  tab <- bind_rows(tab, df)
}

tab_wide <- tab |> 
  pivot_wider(names_from = sample, values_from = Raw_counts)


# substitution of the column names from the .tsv file to the case ID
# re-naming
tab_wide_renamed <- tab_wide |> 
  pivot_longer(cols = contains('.tsv'), names_to = 'File Name', values_to = 'Counts') |> 
  left_join(samplesheet_TCGA_LMS, by = 'File Name') |> 
  dplyr::select(-`File Name`) |> 
  pivot_wider(names_from = 'Case ID', values_from = 'Counts')


# DESEq to have the normalized counts, then to compute the zscore across gene of interest
coldata <- tibble(samples = colnames(tab_wide_renamed |> column_to_rownames("symbol")))
dds <- DESeqDataSetFromMatrix(countData = tab_wide_renamed |> column_to_rownames("symbol"),
                              colData = coldata,
                              design = ~ 1)
dds <- DESeq(dds)
vst <- vst(dds, blind = T)
vst_counts <- as.data.frame(assay(vst)) |> 
  rownames_to_column('symbol')

# save csv
write_csv(vst_counts, "vst_counts_TGCA_leio.csv")


##########
# TCGA ALL
##########

# read samplesheet, we keep only the important infos
samplesheet_TCGA_ALL <- read_delim(samplesheet_TCGA_ALL)

# It is normal for TCGA data to have one patient with multiple RNAseq file
# We filter keeping only primary tumor (i forgot to NOT donwload the metastatic) and that is first aliquot (-01A)
# there are still some duplicates but is normal, so we keep just one of the duplicates
samplesheet_TCGA_ALL <- samplesheet_TCGA_ALL %>% 
  dplyr::filter(`Tumor Descriptor` == 'Primary' & str_detect(`Sample ID`, '-01A')) %>% 
  group_by(`Case ID`) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(`File Name`, `Case ID`)


if (any(duplicated(samplesheet_TCGA_ALL$`Case ID`))) {
  print('There are duplicates Case ID, check the samplesheet')
} else {
  print('Okay, No duplicates :)')
}

# to read, sort and bind in one big tab this is way faster than a for loop
TCGA_all_tab <- map_dfr(files, function(x) {
  read_delim(x,  comment = '#', skip = 6, col_names = F) %>%
    dplyr::select(X2, X4) |> 
    dplyr::rename('symbol' = 'X2', 'Raw_counts' = 'X4') |> 
    mutate(sample = basename(x)) |> 
    dplyr::filter(!duplicated(symbol))
})

TCGA_all_tab_wide <- TCGA_all_tab |> 
  pivot_wider(names_from = sample, values_from = Raw_counts)

# substitution of the column names from the .tsv file to the case ID
# re-naming
TCGA_all_tab_wide_renamed <- TCGA_all_tab_wide |> 
  pivot_longer(cols = contains('.tsv'), names_to = 'File Name', values_to = 'Counts') |> 
  inner_join(samplesheet_TCGA_ALL, by = 'File Name') |> 
  dplyr::select(-`File Name`) |> 
  dplyr::mutate(Counts = ifelse(is.na(Counts), 0, Counts)) %>% 
  pivot_wider(names_from = 'Case ID', values_from = 'Counts') %>% 
  dplyr::filter(!is.na(symbol)) %>% 
  column_to_rownames("symbol")

# DESEq to have the normalized counts, then to compute the zscore across gene of interest
coldata <- tibble(samples = colnames(TCGA_all_tab_wide_renamed))
dds <- DESeqDataSetFromMatrix(countData = TCGA_all_tab_wide_renamed,
                              colData = coldata,
                              design = ~ 1)
dds <- DESeq(dds)
vst <- vst(dds, blind = T)
vst_counts <- as.data.frame(assay(vst)) |> 
  rownames_to_column('symbol')

write_csv(vst_counts, "vst_counts_TGCA_all.csv")
