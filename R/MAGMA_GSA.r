#' Run MAGMA GSA using gene lists of up and downregulated genes obtained using run_diffExp_analysis.r
#'
#' \code{MAGMA_GSA} Uses a dif. expression analysis file produced using run_diffExp_analysis.r to obtain lists of up and down-regulated genes to be used as gene sets and test for overrepresentation in GWAS summary statistics.
#'
#' @param tt_filename diff. Expression file obtained using run_diffExp_analysis function.
#' @param analysis_type Accepted values are 'H_MN-AH', 'D_MN-AH', 'H_Soma-Axon', 'D_Soma-Axon', 'Soma_als-ctrl', or 'Axon_als-ctrl'. Same as indicated in run_diffExp_analysis function.
#' @param genes.raw_path genes.raw path including filename. genes.raw file is obtained using map_snps_togenes formula.
#' @param species Accepted values include "human" or "mouse" default option is mouse.
#' @param gene_n Number of genes to be considered up or downregulated. Gene selection is based on t-value from the differential expression analysis.
#'
#' @return Creates Genecovariate files of gene_n up and downregulated genes, and final GSA of both sets of genes in the given GWAS summary statistics.
#'
#' @examples
#' tt_filename <- run_diffExp_analysis(mouse_D_Soma-Axon_annot,mouse_D_Soma-Axon_data,"/Soma-Axon_Data/","D_Soma-Axon","mouse")
#' MAGMA_GSA(tt_filename = tt_filename, analysis_type = "D_Soma-Axon", genes.raw_path = "/magmafiles/als.genes.raw", species = "mouse", gene_n = 250)
#'
MAGMA_GSA <- function(tt_filename, analysis_type, genes.raw_path, species = "mouse", gene_n){
  if(!analysis_type %in% c("H_MN-AH","D_MN-AH","H_Soma-Axon","D_Soma-Axon","Soma_als-ctrl","Axon_als-ctrl","MN_als-ctrl","AH_als-ctrl")){stop("Incorrect analysis type specified. Analysis_type argument must be set to 'H_MN-AH','D_MN-AH','H_Soma-Axon','D_Soma-Axon','Soma_als-ctrl', or'Axon_als-ctrl'")}
  if(species=="mouse"){
		GetDRList = sprintf("awk -F, '{print $9,$5,$6}' '%s' |sort -g -k 2 |cut -d' ' -f 1 | head -%i > Results_%s_%s/Tables/%s_downreg.txt",tt_filename,gene_n+1,analysis_type,species,gene_n)
		system(GetDRList)
		printf("file %s_downreg.txt stored in the Tables subdirectory of results. \n",gene_n)
		GetURList = sprintf("awk -F, '{print $9,$5,$6}' '%s' |sort -g -k 2 |cut -d' ' -f 1 | sed -e 1b -e :a -e '$q;N;%i,$D;ba' > Results_%s_%s/Tables/%s_upreg.txt",tt_filename,gene_n+2,analysis_type,species,gene_n)
		system(GetURList)
		printf("file %s_upreg.txt stored in the Tables subdirectory of results. \n",gene_n)
	}else{
		GetDRList = sprintf("awk -F, '{print $8,$4,$5}' '%s' |sort -g -k 2 |cut -d' ' -f 1 | head -%i > Results_%s_%s/Tables/%s_downreg.txt",tt_filename,gene_n+1,analysis_type,species,gene_n)
		system(GetDRList)
		printf("file %s_downreg.txt stored in the Tables subdirectory of results. \n",gene_n)
		GetURList = sprintf("awk -F, '{print $8,$4,$5}' '%s' |sort -g -k 2 |cut -d' ' -f 1 | sed -e 1b -e :a -e '$q;N;%i,$D;ba' > Results_%s_%s/Tables/%s_upreg.txt",tt_filename,gene_n+2,analysis_type,species,gene_n)
		system(GetURList)
		printf("file %s_upreg.txt stored in the Tables subdirectory of results. \n",gene_n)
	}
	if(colnames(read.csv(sprintf("Results_%s_%s/Tables/%s_upreg.txt",analysis_type,species,gene_n),sep="\t"))!="HGNC.symbol"){stop("Colnames don't match required format. Check that species indicated is correct")}
	upreg_list <- read.csv(sprintf("Results_%s_%s/Tables/%s_upreg.txt",analysis_type,species,gene_n),sep="\t")[,"HGNC.symbol"]
	downreg_list <- read.csv(sprintf("Results_%s_%s/Tables/%s_downreg.txt",analysis_type,species,gene_n),sep="t")[,"HGNC.symbol"]

	# Write covar files for up and downregulated genes
	write_gene_covarfile(upreg_list,sprintf("UR_%s_%s",analysis_type,species),sprintf("Results_%s_%s/Tables/GeneCovar_UR_%s_%s_%s",analysis_type,species,analysis_type,gene_n,species))
	write_gene_covarfile(downreg_list,sprintf("DR_%s_%s",analysis_type,species),sprintf("Results_%s_%s/Tables/GeneCovar_DR_%s_%s_%s",analysis_type,species,analysis_type,gene_n,species))

	#Map snp to genes and obtain gene-level p-value from summary statistics
	genes_raw_filename <- genes.raw_path

	#GSA using n up-regulated genes
	magma_geneset_UR <- sprintf("magma --gene-results '%s' --set-annot '%s' --out '%s.%s'",genes_raw_filename,sprintf("Results_%s_%s/Tables/GeneCovar_UR_%s_%s_%s",analysis_type,species,analysis_type,gene_n,species),sprintf("Results_%s_%s/Tables/UR_%s_%s",analysis_type,species,analysis_type,gene_n),species)
	system(magma_geneset_UR)

	#GSA using n down-regulated genes
	magma_geneset_DR <- sprintf("magma --gene-results '%s' --set-annot '%s' --out '%s.%s'",genes_raw_filename,sprintf("Results_%s_%s/Tables/GeneCovar_DR_%s_%s_%s",analysis_type,species,analysis_type,gene_n,species),sprintf("Results_%s_%s/Tables/DR_%s_%s",analysis_type,species,analysis_type,gene_n),species)
	system(magma_geneset_DR)
}

