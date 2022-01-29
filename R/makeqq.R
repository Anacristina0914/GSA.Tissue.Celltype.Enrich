makeqq <- function(tt_file,sep=",",n_genes){
  library(ggplot2)
  library(ggrepel)
  library(tools)
  library(tidyverse)
  # Read tt data and rename X column as hgnc_symbol MN-Ah
  tt <- read.csv(file = tt_file,sep=sep)
  tt <- tt %>% dplyr::rename(hgnc_symbol = X)

  #Soma-Axon
  tt <- read.csv(file = tt_file,sep=sep)
  tt_filt <- tt %>% drop_na(baseMean,pvalue)
  tt <- tt_filt %>% dplyr::rename(hgnc_symbol = HGNC.symbol)

  # Read GWAS data (Both Soma-Axon and MN-AH)
  GWAS_ALS <- read.csv("/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Sumstats_Jan/AIBS_M1_MAGMA_Analysis/MAGMA_Files/ALS_sumstats_EUR_formatted.txt.10UP.1.5DOWN/ALS_sumstats_EUR_formatted.txt.10UP.1.5DOWN.genes_whgnc.out",sep="\t")

  # Merge tt data to contain probability of DE(in tt file) v.s. GWAS adjusted P-value (in GWAS file).
  GWAS_ALS <- GWAS_ALS %>% dplyr::select(GWAS_P,hgnc_symbol)
  tt_GWAS_enrich <-  merge(GWAS_ALS,tt,by="hgnc_symbol")

  #For MN-AH analysis
  tt_GWAS_sorted <-  tt_GWAS_enrich[order(tt_GWAS_enrich$GWAS_P,tt_GWAS_enrich$adj.P.Val),]

  # For Soma-Axon
  tt_GWAS_sorted <-  tt_GWAS_enrich[order(tt_GWAS_enrich$GWAS_P,tt_GWAS_enrich$padj),]
  tt_GWAS_sorted$padj_calc <- p.adjust(tt_GWAS_sorted$pvalue, method="fdr")

  # For both Soma-Axon and MN-AH
  plot_source <- file_path_sans_ext(basename(tt_file))
  file_location <- dirname(dirname(tt_file))
  utils::write.table(tt_GWAS_sorted,file=paste(file_location,"/Tables/",plot_source,"GWAS_mergedmod.txt",sep=""),quote=FALSE,row.names=FALSE,sep="\t",col.names=TRUE)
  pdf(paste(file_location,"/Figures/",plot_source,"_DEGWAS1mod.pdf",sep=""),width=11)
  op <- par(mar=c(8,4,4,4))

  # For MN-AH analysis
  tt_enrich_plot <- ggplot(tt_GWAS_sorted[1:n_genes,], aes(x = -log(adj.P.Val), y = -log(GWAS_P))) + geom_point()
  tt_enrich_plot + geom_hline(yintercept = -log(5e-8),linetype="dashed") + geom_vline(xintercept = -log(0.05),linetype="dashed") +
  geom_point(aes(colour = cut(logFC, c(-Inf,-1, 1,Inf))),
               size = 2) +
  scale_color_manual(name ="logFC",
                       values = c("(-Inf,-1]" = "steelblue",
                                  "(-1,1)" = "gray",
                                  "[1, Inf]" = "red"),
                       labels = c("logFC <= -1", "-1 < logFC < 1", "logFC>= 1")) +
  geom_text_repel(aes(label = hgnc_symbol), size = 2) +  labs(title="Adjusted probability DE v.s. GWAS enrichment",caption = paste("Data source: ",plot_source,sep=""))
  dev.off()

  # For Soma-Axon analysis
  tt_enrich_plot <- ggplot(tt_GWAS_sorted[1:n_genes,], aes(x = -log(padj_calc), y = -log(GWAS_P))) + geom_point()
  tt_enrich_plot + geom_hline(yintercept = -log(5e-8),linetype="dashed") + geom_vline(xintercept = -log(0.05),linetype="dashed") +
  geom_point(aes(colour = cut(log2FoldChange, c(-Inf,-1, 1,Inf))),
        size = 2) +
  scale_color_manual(name ="logFC",
                  values = c("(-Inf,-1]" = "steelblue",
                              "(-1,2)" = "gray",
                              "[1, Inf]" = "red"),
                  labels = c("logFC <= -1", "-1 < logFC < 1", "logFC>= 1")) +
  geom_text_repel(aes(label = hgnc_symbol), size = 2) +  labs(title="Adjusted probability DE v.s. GWAS enrichment",caption = paste("Data source: ",plot_source,sep=""))
  dev.off()
}


