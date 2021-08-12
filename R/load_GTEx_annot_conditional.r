#' Load annotation data for MN-AntHorn or Soma-Axon datasets.
#'
#' \code{load_GTEX_annot_conditional} Loads MN-AntHorn or Soma-Axon annotation for matching dataset loaded using the load_GTEX_data_conditional function.
#'
#' @param path path and file name of the GTEx expression dataset. rows=gene_ids, cols=samples.
#' @param data dataset containing expression data, loaded using the load_GTEX_data_conditional function.
#' @param data_type Takes two possible valid options "MN-AH" or "Soma-Axon", to specify which dataset is being loaded, either the Motor Neuron-Anterior horn (GSE18920) or the Soma-Axon (GSE121069) datasets.
#'
#' @return an annotation dataframe containing 3 columns (status: ALS or CTRL, gender: F or M, age: numerical parameter) for the MN-AH, or 2 columns (status: ALS, CTRL, or CTRL_Single; area: Axon, Soma or Single) for the Soma-Axon dataset.
#'
#' @examples
#' # Loads the annotation dataframe for the soma-axon GSE121069_GEO_rpkms_human.txt dataset.
#' annot_human_soma_axon <- load_GTEX_annot_conditional("GSE121069_GEO_rpkms_human.txt",data_human_soma_axon,data_type="Soma-Axon")

load_GTEX_annot_conditional <- function(path,data,data_type="MN-AH"){
	if(!data_type %in% c("MN-AH","Soma-Axon")){stop("This function can only process limited datasets. data_type argument must be set to either 'MN-AH' or 'Soma-Axon'")}
	if(data_type=="MN-AH"){
		annot_spinal_status = rep("ALS",dim(data)[2])
		annot_spinal_status[grep("CTRL",colnames(data))] = "CTRL"
		annot_spinal_area=rep("MN",dim(data)[2])
		annot_spinal_area[grep("AH",colnames(data))] = "AH"
		annot_spinal_gender = rep("Male",dim(data)[2])
		annot_spinal_gender = rep("Male",dim(data)[2])
		annot_spinal_gender[grep("_F",colnames(data))] = "Female"
		annot_spinal_age = as.numeric(as.character(gsub("M|F","",gsub("No.*","",gsub(".*_","",colnames(data))))))
		annot_spinal = data.frame(status=annot_spinal_status,gender=annot_spinal_gender,age=annot_spinal_age,area=annot_spinal_area)
	}else{
		annot_spinal_status = rep("ALS",dim(data)[2])
		annot_spinal_status[grep("Ctrl",colnames(data))] = "CTRL"
		annot_spinal_status[grep("Human",colnames(data))] = "CTRL_Single"
		annot_spinal_area=rep("Axon",dim(data)[2])
		annot_spinal_area[grep("soma",colnames(data))] = "Soma"
		annot_spinal_area[grep("single",colnames(data))] = "Single"
		annot_spinal = data.frame(status=annot_spinal_status,area=annot_spinal_area)
	}
	return(annot_spinal)
}
