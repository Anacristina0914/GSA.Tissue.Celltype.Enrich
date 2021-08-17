# GSA.Tissue.Celltype.Enrich

Ana Cristina Gonzalez Sanchez<br/>
University of Bologna / Karolinska Institutet<br/>
[Contact](mailto:ana.gonzalezsanchez@studio.unibo.it)
 
- [Introduction](#introduction)
- [Installation](#installation)
- [Geneset Enrichment Analysis](#Gene-Set-Enrichment-Analysis-(GSA))
	- [GWAS formatting](#GWAS-Summary-Statistics-Formatting)
	- [Running GSA](#GSA-Enrichment-for-up-and-down-regulated-genes-in-GTEx-data)
- [Tissue Enrichment Analysis](#Tissue-Enrichment-Analysis)
	- [Dataset preparation](#Tissue-Dataset-Preparation)
## Introduction
GSA_Tissue_Celltype_Enrich is an R library intented to carry out in an easy manner Tissue enrichment analysis using data from GTEX, Celltype Enrichment Analysis using MAGMA-Celltyping or Geneset Enrichment Analysis for list of genes coming from differential expression analysis. Over representation of list of genes is tested in GWAS summary statistics. 

## Installation
``` R
if(!"devtools" %in% row.names(installed.packages())){
  install.packages("devtools")
}
library(devtools)

if(!"GSA_Tissue_Celltype_Enrich" %in% row.names(installed.packages())){
  install_github("Anacristina0914/GSA.Tissue.Celltype.Enrich")
}
library(GSA_Tissue_Celltype_Enrich) 
```

## Gene Set Enrichment Analysis (GSA)

### GWAS Summary Statistics Formatting
GWAS summary statistics must be formatted as specified in the [MAGMA](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf) guidance. Briefly, the first three columns have to correspond to the SNP, CHR and BP (see example below):

|SNP |CHR |BP |P |
|:-:|:-:|:-:|:-:|
|rs1345 |2 |100123 |0.01 |
|rs18667 |3 |30566921 |0.5611 |
|rs145 |16 |9992021 |0.173 |

This can be achieved either manually or using the formula format_sumstats_for_magma.r from the [MAGMA_Celltyping](https://github.com/NathanSkene/MAGMA_Celltyping) package.

### GSA Enrichment for up and down regulated genes in GTEx data
``` R
# Load required libraries for the package
load_required_libraries()

# Load GTEx gene expression data from Soma and Axon in human Motor Neurons. Only controls are loaded (C*).  
Human_SomaAxon_data <- GSA.Tissue.Celltype.Enrich::load_GTEx_data_conditional(path = "/Soma_Axon_RNA-Seq/GSE121069_GEO_rpkms_human.txt",pattern = "C*",sep = "\t")

# Load GTEx data annotation
Human_SomaAxon_annot <- GSA.Tissue.Celltype.Enrich::load_GTEx_annot_conditional(path = "/Soma_Axon_RNA-seq/", data = Human_SomaAxon, data_type = "Soma-Axon")

# Run Differential Expression analysis 
tt_human_SomaAxon <- GSA.Tissue.Celltype.Enrich::run_diffExp_analysis(annot = Human_SomaAxon_annot, data = Human_SomaAxon_data, expr_path = "/Soma_Axon_RNA-seq/", analysis_type = "D_Soma-Axon", species = "human")

# Map snps to genes and generate gene level p-value from summary statistics
genes.raw_path <- GSA.Tissue.Celltype.Enrich::map_snps_to_genes(gwas_path = "/ALS_sumstats.txt", N=NULL, genloc_filepath = "/genloc_files/NCBI37.3.gene.loc", genome_ref_path = "/g1000/g1000_eur",analysis_type = "D_Soma-Axon", species = "human")

# Run Gene Set Analysis (GSA) using MAGMA and up/down regulated genes from expression data and GWAS summary statistics.
GSA.Tissue.Celltype.Enrich::MAGMA_GSA(tt_filename = tt_human_SomaAxon, analysis_type = "D_Soma-Axon", genes.raw_path = genes.raw_path, species = "human", gene_n = 250)
```
## Tissue Enrichment Analysis
### Tissue Dataset Preparation
GTEx analysis v8 2017-06-05 gene-level median TPM by tissue dataset was retrieved from the [GTEx database](https://gtexportal.org/home/datasets) and further processed as described in Bryois et al. 2020[[1]](#1) and its corresponding github repository [jbryois/scRNA_disease](https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Code_GTEx/get_GTEx_input.md). The final file used for the tissue enrichment analysis ([top10.txt](https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Code_GTEx/MAGMA/top10.txt)) contains the 10% most specific genes by tissue, where specificity is calculated dividiving the expression of each gene in a tissue divided by total expression of that gene in all tissues.  
``` R
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## References
<a id="1">[1]</a> 
Bryois, J., Skene, N. G., Hansen, T. F., Kogelman, L. J., Watson, H. J., Liu, Z., ... & Sullivan, P. F. (2020). Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson’s disease. Nature genetics, 52(5), 482-493.

## License
[MIT](https://choosealicense.com/licenses/mit/)
