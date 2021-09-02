#' Run differential expression analysis for GTEx data previously loaded using load_GTEx_data_conditional.
#'
#' \code{run_diffExp_analysis} Runs differential expression analysis for GTEx MN-AH or Soma-Axon data.
#'
#' @param annot dataframe containing the annotation dataframe of the expr. data, created using load_GTEx_annot_conditional function.
#' @param data dataframe containing expression data, loaded using the load_GTEX_data_conditional function.
#' @param expr_path Path to the file containing the expression data. Path doesn't include the file name.
#' @param analysis_type string indicating analysis to be carried out. Permitted strings are 'H_AH-MN', 'D_AH-MN', 'H_Axon-Soma', 'D_Axon-Soma', 'Soma_ctrl-als', or 'Axon_ctrl-als'.
#' @param species Permitted strings are either "mouse" or "human".
#'
#' @return path to the differential expression file.
#'
#' @examples
#' # Runs differential expression analysis for mouse Soma-Axon data.
#' tt_mouse_SomaAxon <- run_diffExp_analysis(annot = annot_mouse_axon_soma, data = mouse_axon_soma, expr_path = "/Users/AnaCrisGlez/Soma_Axon_RNA-Seq/", analysis_type = "D_Axon-Soma", species = "mouse")
#'
run_diffExp_analysis <- function(annot, data, expr_path, analysis_type, species = "human"){
	if(!analysis_type %in% c("H_AH-MN","D_AH-MN","H_Axon-Soma","D_Axon-Soma","Soma_crtl-als","Axon_ctrl-als","MN_ctrl-als","AH_ctrl-als")){stop("Incorrect analysis type specified. Analysis_type argument must be set to 'H_AH-MN','D_AH-MN','H_Axon-Soma','D_Axon-Soma','Soma_ctrl-als', or'Axon_ctrl-als'")}
	if(!species %in% c("human","mouse")){stop("Incorrect species specified. Species must be set to either mouse or human")}
	setwd(expr_path)
	if (!file.exists(paste("Results_",analysis_type,"_",species,sep = ""))){
		dir.create(paste("Results_",analysis_type,"_",species,sep = ""))
		dir.create(paste("Results_",analysis_type,"_",species,"/Tables",sep = ""))
		dir.create(paste("Results_",analysis_type,"_",species,"/Figures",sep = ""))
	}
	file_name = paste("tt_",analysis_type,".csv",sep="")
	if(analysis_type == "H_Axon-Soma" || analysis_type == "D_Axon-Soma"){
		mod  = model.matrix(~annot$area,levels = c("Axon","Soma"))
		colnames(mod)[2] = "Area"
		coef = "Area"
	} else if (analysis_type == "Axon_ctrl-als" || analysis_type == "Soma_ctrl-als"){
		mod  = model.matrix(~annot$status,levels = c("CTRL","ALS"))
		colnames(mod)[2] = "Status"
		coef = "Status"
	} else if (analysis_type == "MN_ctrl-als" || analysis_type == "AH_ctrl-als"){
	  mod  = model.matrix(~annot$gender+factor(annot$status,levels=c("CTRL","ALS")))
	  colnames(mod)[2:3] = c("Gender","Status")
	  print(mod)
	  coef = "Status"
	} else if ((analysis_type == "H_AH-MN" || analysis_type == "D_AH-MN") && species == "human"){
		mod  = model.matrix(~annot$gender+factor(annot$area,levels=c("AH","MN")))
		colnames(mod)[2:3] = c("Gender","Area")
		coef = "Area"
	} else{
		mod = model.matrix(~annot$area, levels = c("AH","MN"))
		colnames(mod)[2] = "Area"
		coef = "Area"
	}
tt_spinal = prep.tt(data,mod,coef = coef,species = species,analysis_type = analysis_type,file_name = file_name)
  return(paste(getwd(),sprintf("/Results_%s_%s/Tables/%s",analysis_type,species,file_name),sep=""))
}
