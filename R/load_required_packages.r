#' Load required packages to run GSA_Tissue_Celltype_Enrichment.
#'
#' \code{load_required_libraries} Load required packages to run GSA_Tissue_Celltype_Enrichment, if packages are not installed, the code installs them before loading them.
#'
#' @examples
#' load_required_libraries()

load_required_libraries <- function(){
	required_packages_library=c("devtools","R.utils","tidyr","limma","RNOmni","dplyr","tidyverse")
	required_packages_github=c("MAGMA.Celltyping","One2One","EWCE")
	missing_index_library <- which(!required_packages_library %in% row.names(installed.packages()),TRUE)
	missing_index_github <- which(!required_packages_github %in% row.names(installed.packages()),TRUE)
	for(i in missing_index_library){
		install.packages(required_packages[i])
	}
	library(devtools)
	for (j in missing_index_github){
		if(required_packages_github[j]=="MAGMA.Celltyping"){
		devtools::install_github("NathanSkene/MAGMA_Celltyping")
		} else if(required_packages_github[j]=="One2One"){
		devtools::install_github("NathanSkene/One2One")
		} else{
		devtools::install_github("neurogenomics/EWCE")
		}
	}
	library(MAGMA.Celltyping)
	library(One2One)
	library(R.utils)
	library(tidyr)
	library(limma)
	library(EWCE)
	library(RNOmni)
	library(dplyr)
	library(tidyverse)
}
