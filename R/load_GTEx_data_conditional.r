#' Load GTEx MN-AntHorn or Soma-Axon datasets.
#'
#' \code{load_GTEX_data_conditional} Loads MN-AntHorn or Soma-Axon datasets from GTEx in the format rows=gene_ids, cols=samples. Assumes that the first element of the dataset is "NAME".
#'
#' @param path path and file name of the GTEx expression dataset. rows=gene_ids, cols=samples.
#' @param pattern grep patter to match samples to be kept for the analysis.
#' @param sep sep argument to denote separation used in the expression dataset. Valid options are comma="," tab="\t" space=" " semicolon=";" etc.
#'
#' @return a dataframe containing n rows equal to genes and m columns equal to the samples matching the grep pattern specified.
#'
#' @examples
#' # Loads the control samples of soma and axon within the GSE121069_GEO_rpkms_human.txt dataset.
#' data_human_soma_axon <- load_GTEX_data_condition("GSE121069_GEO_rpkms_human.txt","C*",sep="\t")
#' grep pattern_example (for file GSE18920): "*_ALS_*" = loads only cases (MN and AH), *_CTRL_*" = loads only controls (MN and AH).
load_GTEX_data_conditional <- function(path,pattern,sep="\t"){
	data_spinal = read.csv(path,stringsAsFactors = FALSE,sep=sep)
	rownames(data_spinal) = gsub(" ","",data_spinal$NAME)
	#Eliminates spaces in gene_ids in rows.
	data_spinal = data_spinal[,-1]
	data_spinal = data_spinal[,grep(pattern,colnames(data_spinal))]
	#Keeps columns that match the given grep sample pattern.
	return(data_spinal)
}
