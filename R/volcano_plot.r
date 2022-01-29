if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library("EnhancedVolcano")
library("tidyverse")

#Change Dif_expfile name to generate the new plot
# All Soma-Axon DEA
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Soma_Axon_RNA-Seq/Results_Soma_als-ctrl_mouse/Tables/tt_Soma_als-ctrl.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Soma_Axon_RNA-Seq/Results_H_Soma-Axon_mouse/Tables/tt_H_Soma-Axon.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Soma_Axon_RNA-Seq/Results_H_Soma-Axon_human/Tables/tt_H_Soma-Axon.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Soma_Axon_RNA-Seq/Results_Axon_als-ctrl_mouse/Tables/tt_Axon_als-ctrl.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/Soma_Axon_RNA-Seq/Results_D_Soma-Axon_mouse/Tables/tt_D_Soma-Axon.csv"

# MN-AH DEA
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/MotNeur-AntHorn_Data/Results_AH_als-ctrl_human/Tables/tt_AH_als-ctrl.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/MotNeur-AntHorn_Data/Results_D_MN-AH_human/Tables/tt_MN-AH.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/MotNeur-AntHorn_Data/Results_H_MN-AH_human/Tables/tt_MN-AH.csv"
Dif_expfile <- "/Users/AnaCrisGlez/Desktop/Bioinformatics MSc./Erasmus+/Thesis/Data/MotNeur-AntHorn_Data/Results_MN_als-ctrl_human/Tables/tt_MN_als-ctrl.csv"

tt_file <- read.csv(Dif_expfile)
file_location <- dirname(dirname(Dif_expfile))

if(grepl("human",file_location)==TRUE){org <- "human"}else{org <- "mouse"}
plot_name <- paste(file_path_sans_ext(basename(Dif_expfile)),"_",org,sep="")

if(grepl("Axon",plot_name)==TRUE||grepl("Soma",plot_name)==TRUE){
  log2 <- "log2FoldChange"
  pval <- "pvalue"
}else{
  log2 <- "logFC"
  pval <- "P.Value"
}

pdf(paste(file_location,"/Figures/",plot_name,".pdf",sep=""),width=11)
op <- par(mar=c(8,4,4,4))
EnhancedVolcano(tt_file,lab=tt_file$HGNC.symbol,x=log2,y=pval,title=plot_name,shape = c(1, 4, 23, 25),colAlpha = 1,labSize = 4,drawConnectors = TRUE,widthConnectors = 0.4)
dev.off()

