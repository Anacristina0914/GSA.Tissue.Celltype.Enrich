#' Convert an expression dataset in rpkm to tpm. 
#' 
#' \code{rpkm_to_tpm} convert a dataset expressed in rpkm into tpm. 
#' @param dataset expressed in rpkm. rows = genes and columns = samples. 
#' 
#' @return a dataset expressed in tpm.
#'
#' @example 
#' data_human_soma_axon <- load_GTEX_data_conditional("GSE121069_GEO_rpkms_human.txt","C*",sep="\t")
#' data_human_rpkm <- rpkm_to_tpm(data_human_soma_axon)

rpkm_to_tpm <- function(dataset){
  TPM <- apply(dataset,2,function(x){x/sum(x)*1000000}) 
  return(TPM)
}
