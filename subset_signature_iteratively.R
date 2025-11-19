{
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(Hmisc)
  library(patchwork)
  library(gridExtra)
}

plot_dir <- '/media/alessio/Data/Public_data_TCGA_ENCODE_cBioportal/plot_cluster_survival/'

vst_leio <- 'vst_counts_TGCA_leio.csv'
clinical_leio <- 'manual_downloaded_samples/clinical.tsv'
Signature_33 <- read_lines('/media/alessio/Data/hypoxia/DBR_analysis_RNAseq_new/Genes_Up_shared_DBR_promoter.txt')

vst_file <- vst_leio
clinical_file <- clinical_leio
genes <- Signature_33


# prepare data
{
  ### read and compute Zscore
  vst_counts <- read_delim(vst_file) |>
    column_to_rownames('symbol')
  
  # Zscore
  Zscore_tab <- t(scale(t(vst_counts)))
  
  genes <- genes[genes %in% row.names(Zscore_tab)]
  # select genes

}

{
  iterations <- 1
  n_genes <- 4
  centers <- 2
  best_data <- list()
  
  for (i in 1:iterations) {
    
    genes_subset <- sample(genes, n_genes)
    genes_subset <- c('PER2', 'TEF', 'SP110', 'PDE4D')
    namex <- paste(genes_subset, collapse = '_')
    
    Zscore_tab_subset <- Zscore_tab[genes_subset, ] %>% 
      na.omit() %>%  
      as.data.frame() 
    
    Zscore_HDAC4 <-  Zscore_tab['HDAC4', , drop = F] %>%  
      as.data.frame()
    
    ### Kmeans
    # first transpose because we want to clusterize samples
    Zscore_transposed <- t(Zscore_tab_subset)
    
    # result is not clear bu I arbitrary choose 3 :)
    kmeans_samples <- kmeans(Zscore_transposed, centers = centers, nstart = 25)
    
    # Ensure the sample order in the data matches the order in the annotation
    ordered_samples <- rownames(Zscore_transposed)[order(kmeans_samples$cluster)]
    
    # create the annotation
    annotation_row <- data.frame(cluster_samples = factor(kmeans_samples$cluster[ordered_samples]))
    rownames(annotation_row) <- ordered_samples
    
    # Re order the matrix
    Zscore_transposed_ordered <- Zscore_transposed[ordered_samples, ]
    
    # Heatmap
    
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
    
    Zscore_HDAC4 <- t(Zscore_HDAC4)
    Zscore_HDAC4 <- Zscore_HDAC4[ordered_samples, , drop = F]
    
    Zscore_transposed_ordered_with_HDAC4 <- cbind(Zscore_transposed_ordered, Zscore_HDAC4)
    
    
    H <- pheatmap(Zscore_transposed_ordered_with_HDAC4,
                  silent = T,
                  show_colnames = TRUE,
                  show_rownames = FALSE,
                  gaps_col = ncol(Zscore_transposed_ordered), 
                  cluster_rows = F,    
                  cluster_cols = F,   
                  annotation_row = annotation_row,
                  border_color = 'grey',
                  annotation_colors = annotation_colors,
                  main = 'Zscore Expression by K-means Cluster',
                  color = myColor,
                  breaks = myBreaks,
                  cellwidth = 10
    )
  
    
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
    res <- survdiff(Surv(time, status) ~  cluster_samples, data = clinical_unique)
    
    surv <- ggsurvplot(
      fit,
      pval = TRUE,
      conf.int = FALSE,
      risk.table = FALSE,
      legend.title = '',
      palette = c("blue", "red"),
      ggtheme = theme_bw(base_size = 15) + 
        theme(aspect.ratio = 0.5, 
              plot.margin = margin(10, 10, 10, 10))
    )
    surv_grob <- ggplotGrob(surv$plot)
    
    # extraction of the best cluster
    # take pvalue
    p.val <- 1 - pchisq(res$chisq, length(res$n) - 1)
    
    # take the cluster with better survival
    data <- summary(fit)$table %>% as.data.frame()
    best_cluster_survival <- data[which.max(data$median), ] %>% rownames()
    
    # check if it corresponds with the data with the highest expression
    Zscore_transposed_ordered_1 <- Zscore_transposed_ordered[annotation_row$cluster_samples == 1, ] 
    Zscore_transposed_ordered_2 <- Zscore_transposed_ordered[annotation_row$cluster_samples == 2, ] 
    #Zscore_transposed_ordered_3 <- Zscore_transposed_ordered[annotation_row$cluster_samples == 3, ]
    
    
    # Now we want to keep the result if:
    # - the minimum number of the element in the cluster is 25 for the calculation
    # - the cluster with best survival is also the one with higher expression
    # - the pvalue is less than 0.05
    
    min_number <- 25
    
    means <- tibble(
      cluster = c('cluster_samples=1',
                  "cluster_samples=2"),
      mean = c(
        mean(Zscore_transposed_ordered_1),
        mean(Zscore_transposed_ordered_2)
      ),
      nsamples = c(
        nrow(Zscore_transposed_ordered_1),
        nrow(Zscore_transposed_ordered_2)
      )
    )
    
    max_expression_cluster <- means[which.max(means$mean), ]$cluster
    max_expression_mean <-  means[which.max(means$mean), ]$mean
    
    if (best_cluster_survival == max_expression_cluster & 
        p.val < 0.0001 & 
        !any(means$nsamples < 25)
        ) {
      print('Cluster with max expression has also the best survival \n')
      best_data_l <- list(
        genes = genes_subset,
        cluster = best_cluster_survival,
        mean_expr = max_expression_mean,
        pvalue = p.val
      )
      best_data[[i]] <- best_data_l
      
      grid.arrange(H$gtable, surv_grob, ncol = 2, heights = c(4, 0.5))
      
      pdf(file = paste0(plot_dir,'/',namex,'.pdf'), height = 8, width = 12)
      grid.arrange(H$gtable, surv_grob, ncol = 2, heights = c(4, 0.5))
      dev.off()
      
    } else {
      print('NO correspondence\n Skip \n')
    }
    
  }
}







