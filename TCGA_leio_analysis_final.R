{
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(Hmisc)
}


vst_leio <- 'vst_counts_TGCA_leio.csv'
vst_all <-  'vst_counts_TGCA_all.csv'

clinical_leio <- 'manual_downloaded_samples/clinical.tsv'
clinical_all <- 'all_TCGA/clinical.tsv'

cht_leio <- 5
cht_all <- 0.2

HDACs_MEFs <- c('MEF2B', 'HDAC7', 'MEF2C', 'HDAC4', 'MEF2A', 'HDAC5', 'MEF2D', 'HDAC9')
Signature_33 <- read_lines('/media/alessio/Data/hypoxia/DBR_analysis_RNAseq_new/Genes_Up_shared_DBR_promoter.txt')
Signature_33_plus <- c(read_lines('/media/alessio/Data/hypoxia/DBR_analysis_RNAseq_new/Genes_Up_shared_DBR_promoter.txt'), 'NFIA', 'NFIC', 'NFIX')

vst_file <- vst_leio
clinical_file <- clinical_leio
cht <- cht_leio
genes <- Signature_33
genes_2 <- HDACs_MEFs
heatmap_zscore_name <- 'Leio_Hmap_signature_33.pdf'
survival_plot_name <- 'Leio_Survival_curve_33_genes.pdf'

# Analysis over LMS only TCGA samples
## Compute Zscore, clustering with kmeans and Heatmap

### read and compute Zscore
vst_counts <- read_delim(vst_file) |>
  column_to_rownames('symbol')

# Zscore
Zscore_tab <- t(scale(t(vst_counts)))


genes <- genes[genes %in% row.names(Zscore_tab)]
# select genes
Zscore_tab_subset <- Zscore_tab[genes, ] %>% 
  na.omit() %>%  
  as.data.frame() 

### Kmeans
# first transpose because we want to clusterize samples
Zscore_transposed <- t(Zscore_tab_subset)

# elbow
wss <- numeric()
k_values <- 2:10

for (k in k_values) {
  set.seed(10)
  km <- kmeans(Zscore_transposed,
               centers = k,
               nstart = 25)
  wss[k-1] <- km$tot.withinss
}

plot(k_values, wss, type = 'b', pch = 19)

# result is not clear bu I arbitrary choose 3 :)
kmeans_samples <- kmeans(Zscore_transposed, centers = 2, nstart = 25)

# clusters of samples 
kmeans_samples$cluster

# Ensure the sample order in the data matches the order in the annotation
ordered_samples <- rownames(Zscore_transposed)[order(kmeans_samples$cluster)]

# create the annotation
annotation_row <- data.frame(cluster_samples = factor(kmeans_samples$cluster[ordered_samples]))
rownames(annotation_row) <- ordered_samples

# Re order the matrix
Zscore_transposed_ordered <- Zscore_transposed[ordered_samples, ]


# Heatmap
# Compute distance and hierarchical clustering for columns (genes)
dist_cols <- dist(t(Zscore_transposed_ordered))
hc_cols <- hclust(dist_cols, method = "complete")

# Cut into 3 clusters
colclust <- cutree(hc_cols, k = 3)
# Create annotation dataframe
annotation_col <- data.frame(cluster_genes = factor(colclust))

# Reorder columns by cluster
ordered_cols <- names(sort(colclust))
Zscore_transposed_ordered <- Zscore_transposed_ordered[, ordered_cols]

# Update annotation order to match
annotation_col <- annotation_col[ordered_cols, , drop = FALSE]

# Reading clinical info
clinical <- read_tsv(clinical_file)

# select only the variables of interests
clinical_unique <- clinical |> 
  dplyr::filter(!duplicated(cases.case_id)) |> 
  dplyr::select(cases.submitter_id, demographic.days_to_death, diagnoses.days_to_last_follow_up, demographic.vital_status, cases.primary_site)

# add tissue of origin info 
annotation_row <- annotation_row %>% 
  rownames_to_column('cases.submitter_id') %>% 
  left_join(clinical_unique %>% 
              dplyr::select(cases.submitter_id, cases.primary_site), by = 'cases.submitter_id') %>% 
  dplyr::select(cases.submitter_id, cluster_samples, cases.primary_site) %>% 
  column_to_rownames('cases.submitter_id')

annotation_row$cases.primary_site <- factor(annotation_row$cases.primary_site)

# Plot the heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(4/paletteLength, 4, length.out=floor(paletteLength/2)))

annotation_colors <- list(
  cluster_samples = c('1' = "blue", '2' ="red"),
  cluster_genes = c('1' = 'brown', '2' = 'darkgreen', '3' = 'purple'),
  cases.primary_site = c(
    'Bones, joints and articular cartilage of limbs' = '#FF1744',  # neon red/pink
    'Connective, subcutaneous and other soft tissues' = '#00E5FF', # electric cyan
    'Corpus uteri' = '#76FF03',  # neon lime
    'Other and unspecified parts of tongue' = '#D500F9',  # vivid magenta
    'Ovary' = '#FFAB00',  # bright amber / gold
    'Retroperitoneum and peritoneum' = '#1DE9B6',  # aqua green
    'Stomach' = '#FF6D00',  # hot orange
    'Uterus, NOS' = '#304FFE'  # vivid blue
  )
)

pdf(heatmap_zscore_name, width = length(genes) / 1.6, height = nrow(Zscore_transposed_ordered) / 3.5)
pheatmap(Zscore_transposed_ordered,
         show_colnames = TRUE,
         show_rownames = FALSE,
         cluster_rows = F,    
         cluster_cols = T,   
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         border_color = 'grey',
         annotation_colors = annotation_colors,
         main = 'Zscore Expression by K-means Cluster',
         color = myColor,
         breaks = myBreaks,
         cellwidth = 10,
         cellheight = cht
)
dev.off()


# to have the cluster visualization
annotation_row_df <- annotation_row |> 
  as.data.frame() |> 
  rownames_to_column('cases.submitter_id')

# prepare to plot the survival KM
# convert survival time e status
clinical_unique$time <- ifelse(clinical_unique$demographic.days_to_death == '\'--', 
                               clinical_unique$diagnoses.days_to_last_follow_up,
                               clinical_unique$demographic.days_to_death)

clinical_unique$status <- ifelse(clinical_unique$demographic.vital_status == "Dead", 1, 0)

# add the Kmeans group and filter those with no infos of the time
clinical_unique <- clinical_unique |> 
  left_join(annotation_row_df, by = 'cases.submitter_id') |> 
  dplyr::filter(time != '\'--') |> 
  dplyr::mutate(time = as.numeric(time))

clinical_unique$cluster_samples <- factor(clinical_unique$cluster_samples)

# Kaplan meier
fit <- survfit(Surv(time, status) ~ cluster_samples, data = clinical_unique)

p <- ggsurvplot(fit,
           pval = TRUE,
           conf.int = F,
           risk.table = F,
           legend.title = '',
           palette = c("blue", "red", "#5F99EA"),
           ggtheme =   theme_bw(), 
           newpage = F
)
p
ggsave(p$plot, device = 'pdf', filename = survival_plot_name, height = 6, width = 8)

# Relationship between two list of genes
### heatmap of correlation
core_genes <- genes_2
core_genes <- core_genes[core_genes %in% rownames(Zscore_tab)]


total_genes <- c(core_genes, genes) 

# Heatmap correlations between zscore
Zscore_crossed_genes<- Zscore_tab[total_genes, ]

res <- rcorr(t(Zscore_crossed_genes), type = 'spearman')
cor_mat <- res$r
p_mat <- res$P

# filtering out ns pvalue
# cor_mat[p_mat >= 0.05] <- 0

cor_mat_selected <- cor_mat[core_genes, genes]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

H <- pheatmap(cor_mat_selected,
              show_colnames = T,          
              cluster_rows = T,
              cluster_cols = T,
              main = 'Spearman corr',
              color=myColor,
              breaks=myBreaks)

H


