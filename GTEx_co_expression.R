# GTEx co-expression
{
  library(tidyverse)  
  library(Hmisc)
  library(ComplexHeatmap)
}

genes_signature <- read_lines('/media/alessio/Data/hypoxia/DBR_analysis_RNAseq_new/Genes_Up_shared_DBR_promoter.txt')
HDAC_MEFs <- genes_HDAC_MEF <- c('HDAC4', 'HDAC5', 'HDAC7', 'HDAC9', 'MEF2A', 'MEF2B', 'MEF2C', 'MEF2D')

GTEx_TPM <- read_delim('GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct', comment = '#', skip = 2)

GTEx_TPM_sub <- GTEx_TPM %>% 
  dplyr::filter(Description %in% c(genes_signature, HDAC_MEFs)) %>% 
  dplyr::select(-Name) %>% 
  dplyr::rename('SYMBOL' = 'Description') %>% 
  column_to_rownames('SYMBOL') %>% 
  as.matrix()

genes_signature <- genes_signature[genes_signature %in% rownames(GTEx_TPM_sub)]

type <- 'spearman'

res_corr <- rcorr(t(GTEx_TPM_sub), type = type)
res_corr_r <- res_corr$r
res_corr_p <- res_corr$P

res_corr_r[res_corr_p >= 0.05] <- NA

res_cor_sel <- res_corr_r[genes_signature, HDAC_MEFs]

colfunc <- circlize::colorRamp2(c(-1, 0, 1),
                                c("blue", "white", "red"))
Heatmap(res_cor_sel, 
        name = type, 
        rect_gp = gpar(col = "white", lwd = 2),
        col = colfunc, 
        cluster_rows = F, 
        cluster_columns = F,
        column_title = 'GTEx manual ops')


# Spearman corr from coGTEx (https://victortrevino.bioinformatics.mx:8181/cogtex/geneHome.jsp?Version=Base)
Co_GTEx_33_genes_HDAC_MEFs <- read_delim('CoGTEx_Pearson_33_genes_HDACs_MEFs.txt') %>% 
  column_to_rownames('...1') %>% 
  as.matrix()

genes_signature <- genes_signature[genes_signature %in% rownames(Co_GTEx_33_genes_HDAC_MEFs)]

Co_GTEx_33_genes_HDAC_MEFs_sub <- Co_GTEx_33_genes_HDAC_MEFs[genes_signature, HDAC_MEFs]

colfunc <- circlize::colorRamp2(c(-0.5, 0, 0.5),
                                c("blue", "white", "red"))

pdf('CoGTEx_signature_genes_HDAC_MEFs.pdf', width = 5, height = 8)
Heatmap(Co_GTEx_33_genes_HDAC_MEFs_sub, 
        name = type, 
        rect_gp = gpar(col = "white", lwd = 2),
        col = colfunc,
        cluster_rows = F, 
        cluster_columns = F, 
        column_title = 'CoGTEx')
dev.off()
