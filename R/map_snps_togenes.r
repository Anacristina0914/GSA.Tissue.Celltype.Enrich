#' Write .genes.annot and .genes.raw files required for GSA analysis using MAGMA.
#'
#' \code{map_snps_togenes} This function precedes MAGMA_GSA function if now .genes.raw file is present. It is used to prepare GWAS file for GSA.
#'
#' @param gwas_path GWAS summary statistics path including file name.
#' @param N N parameter of GWAS summary statistics can be set to the total N. If no N column is present in GWAS summary statistics, this parameter has to be specified by the user.
#' @param genloc_filepath Human gene location file path.
#' @param genome_ref_path Human genome reference path and prefix of .bed, .bim, .fam files. If .synonyms file with the same prefix is present in the same folder it is automatically detected.
#' @param analysis_type Accepted values are 'H_MN-AH', 'D_MN-AH', 'H_Soma-Axon', 'D_Soma-Axon', 'Soma_als-ctrl', or 'Axon_als-ctrl'. Same as indicated in run_diffExp_analysis function.
#' @param species Gene expression species. Accepted values include "human" or "mouse" default option is mouse.
#'
#' @return Path of the .genes.raw file to be used for MAGMA_GSA. Two files (.genes.annot and .genes.raw) are generated and stored. Maps snps in GWAS summary statistics into genes and produces a gene-level p-value to be used for the GSA.
#'
#' @examples
#' genes.raw_path <- map_snps_togenes(gwas_path = "/Sumstats_ALS/ALS_sumstats_EUR_only.txt",N = NULL,genloc_filepath = "/genloc_files/NCBI37.3.gene.loc", genome_ref_path = "/g1000/g1000_eur", analysis_type = "D_Soma-Axon", species = "mouse")
#'
map_snps_togenes  <- function(gwas_path,N=NULL,genloc_filepath,genome_ref_path,analysis_type,species){
  if(!analysis_type %in% c("H_MN-AH","D_MN-AH","H_Soma-Axon","D_Soma-Axon","Soma_als-ctrl","Axon_als-ctrl","MN_als-ctrl","AH_als-ctrl")){stop("Incorrect analysis type specified. Analysis_type argument must be set to 'H_MN-AH','D_MN-AH','H_Soma-Axon','D_Soma-Axon','Soma_als-ctrl', or'Axon_als-ctrl'")}
  setwd(paste("Results_",analysis_type,species,sep=""))
	outfile=paste(analysis_type,"_",species,sep="")
	snp_to_gene <- sprintf("magma --annotate window=10,1.5 --snp-loc '%s' --gene-loc '%s' --out '%s'",gwas_path,genloc_filepath,outfile)
	system(snp_to_gene)
	if(is.null(N)){
		con <- file(gwas_path,"r") ; first_line <- readLines(con,n=1) ; close(con)
		column_headers = strsplit(first_line,"\t")[[1]]
		if("N" %in% column_headers){n_arg = "ncol=N"}else{
			nval <- as.numeric(readline("There is no N column within the sumstats file. What is the N value for this GWAS?"))
			if(is.na(nval)){stop(sprintf("%s provided but value of N for the GWAS must be numeric",nval))}
			if(nval<1000){stop("Value of N provided is less than 1000. This seems unlikely.")}
			if(nval>100000000){stop("Value of N provided is over than 100000000. In 2018 this seems unlikely.")}
			n_arg = sprintf("N=%s",nval)
		}
	}else{
		n_arg = sprintf("N=%s",N)
	}
	genepval <- sprintf("magma --bfile '%s' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'", genome_ref_path,gwas_path,n_arg,outfile,outfile)
	system(genepval)
	return(paste(getwd(),"/",outfile,".genes.raw",sep=""))
}
