# GSA.Tissue.Celltype.Enrich

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

## Usage

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
```
GWAS summary statistics must be formatted as specified in the [MAGMA](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf) guidance. Briefly, the first three columns have to correspond to the SNP, CHR and BP (see example below):

|SNP |CHR |BP |P |
---------------------------------
|rs1345 |2 |100123 |0.01 |
---------------------------------
|rs18667 |3 |30566921 |0.5611 |
------------------------------------------------
|rs145 |16 |9992021 |0.173 |
----------------------------------------
This can be achieved either manually or using the formula format_sumstats_for_magma.r from the [MAGMA_Celltyping](https://github.com/NathanSkene/MAGMA_Celltyping) package.

``` R
# Map snps to genes and generate gene level p-value from summary statistics
genes.raw_path <- GSA.Tissue.Celltype.Enrich::map_snps_to_genes(gwas_path = "/ALS_sumstats.txt", N=NULL, genloc_filepath = "/genloc_files/NCBI37.3.gene.loc", genome_ref_path = "/g1000/g1000_eur",analysis_type = "D_Soma-Axon", species = "human")

# Run Gene Set Analysis (GSA) using MAGMA and up/down regulated genes from expression data and GWAS summary statistics.
GSA.Tissue.Celltype.Enrich::MAGMA_GSA(tt_filename = tt_human_SomaAxon, analysis_type = "D_Soma-Axon", genes.raw_path = genes.raw_path, species = "human", gene_n = 250)
```
### Tissue Enrichment in GWAS summary statistics
``` R
[GTEx_data_v8]("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz") was retrieved.

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
