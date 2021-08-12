#' Run differential expression analysis for GTEx data previously loaded using load_GTEx_data_conditional.
#'
#' \code{run_diffExp_analysis} Runs differential expression analysis for GTEx MN-AH or Soma-Axon data.
#'
#' @param annot dataframe containing the annotation dataframe of the expr. data, created using load_GTEx_annot_conditional function.
#' @param data dataframe containing expression data, loaded using the load_GTEX_data_conditional function.
#' @param expr_path Path to the file containing the expression data. Path doesn't include the file name.
#' @param analysis_type string indicating analysis to be carried out. Permitted strings are 'H_MN-AH', 'D_MN-AH', 'H_Soma-Axon', 'D_Soma-Axon', 'Soma_als-ctrl', or 'Axon_als-ctrl'.
#' @param species Permitted strings are either "mouse" or "human".
#'
#' @return path to the differential expression file.
#'
#' @examples
#' # Runs differential expression analysis for mouse Soma-Axon data.
#' tt_mouse_SomaAxon <- run_diffExp_analysis(annot = annot_mouse_axon_soma, data = mouse_axon_soma, expr_path = "/Users/AnaCrisGlez/Soma_Axon_RNA-Seq/", analysis_type = "D_Soma-Axon", species = "mouse")
#'
run_diffExp_analysis <- function(annot, data, expr_path, analysis_type, species = "human"){
	if(!analysis_type %in% c("H_MN-AH","D_MN-AH","H_Soma-Axon","D_Soma-Axon","Soma_als-ctrl","Axon_als-ctrl","MN_als-ctrl","AH_als-ctrl")){stop("Incorrect analysis type specified. Analysis_type argument must be set to 'H_MN-AH','D_MN-AH','H_Soma-Axon','D_Soma-Axon','Soma_als-ctrl', or'Axon_als-ctrl'")}
	if(!species %in% c("human","mouse")){stop("Incorrect species specified. Species must be set to either mouse or human")}
	setwd(expr_path)
	if (!file.exists(paste("Results_",analysis_type,"_",species,sep = ""))){
		dir.create(paste("Results_",analysis_type,"_",species,sep = ""))
		dir.create(paste("Results_",analysis_type,"_",species,"/Tables",sep = ""))
		dir.create(paste("Results_",analysis_type,"_",species,"/Figures",sep = ""))
	}
	file_name = paste("tt_",analysis_type,".csv",sep="")
	if(analysis_type == "H_Soma-Axon" || analysis_type == "D_Soma-Axon"){
		mod  = model.matrix(~annot$area,levels = c("Soma","Axon"))
		colnames(mod)[2] = "Area"
		coef = "Area"
	} else if (analysis_type == "Axon_als-ctrl" || analysis_type == "Soma_als-ctrl"){
		mod  = model.matrix(~annot$status,levels = c("ALS","CTRL"))
		colnames(mod)[2] = "Status"
		coef = "Status"
	} else if (analysis_type == "MN_als-ctrl" || analysis_type == "AH_als-ctrl"){
	  mod  = model.matrix(~annot$gender+factor(annot$status,levels=c("ALS","CTRL")))
	  colnames(mod)[2:3] = c("Gender","Status")
	  print(mod)
	  coef = "Status"
	} else if ((analysis_type == "H_MN-AH" || analysis_type == "D_MN-AH") && species == "human"){
		mod  = model.matrix(~annot$gender+factor(annot$area,levels=c("MN","AH")))
		colnames(mod)[2:3] = c("Gender","Area")
		coef = "Area"
	} else{
		mod = model.matrix(~annot$area, levels = c("MN","AH"))
		colnames(mod)[2] = "Area"
		coef = "Area"
	}
tt_spinal = prep.tt(data,mod,coef = coef,species = species,analysis_type = analysis_type,file_name = file_name)
  return(paste(getwd(),sprintf("/Results_%s_%s/Tables/%s",analysis_type,species,file_name),sep=""))
}

