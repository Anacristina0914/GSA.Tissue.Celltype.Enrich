#' Prepare variables and dataset for diff. expression analysis. This function is used by run_als_diffExp_analysis.
#'
#' \code{prep.tt} Function used by run_diffExp_analysis function to prepare data and variables for diff expression analysis of MN-AH or Soma-Axon data.
#'
#' @param exp_data dataframe containing expression data where rows=gene_ids and cols=sample_ids.
#' @annot annot_data dataframe containing annotation data.
#' @param mod Design or model matrix for Diff. Exp. Analysis.
#' @param coef parameter taken by topTable function to indicate which variable of interest is to be extracted from the linear model.
#' @param species species from which expression data was derived. Valid options are either "mouse" or "human".
#' @param analysis_type string indicating analysis to be carried out. Allowed strings are 'H_MN-AH', 'D_MN-AH', 'H_Soma-Axon', 'D_Soma-Axon', 'Soma_als-ctrl', or 'Axon_als-ctrl'.
#' @param file_name Diff. Expression Analysis file name.
#'
#' @return a file is stored in the subfolder Results/Tables of the working directory and a dataframe with the restuls of the diff. expression analysis is returned.

prep.tt <- function(exp_data,annot,mod,coef,species,analysis_type,file_name){
	if(analysis_type %in% c("H_MN-AH","D_MN-AH","MN_als-ctrl","AH_als-ctrl")){
	  fit = lmFit(exp_data,mod)
	  eb = eBayes(fit)
	  tt = topTable(eb, coef=coef, adjust="BH",number=1000000)
	} else{
	    if(analysis_type %in% c("H_Soma-Axon","D_Soma-Axon")){
	      annot$area <- factor(annot$area,levels=c("Soma","Axon"))
	      DE <- DESeqDataSetFromMatrix(countData = exp_data, colData = annot, design= ~area)
	    }else{
	      annot$status <- factor(annot$status,levels=c("ALS","CTRL"))
	      DE <- DESeqDataSetFromMatrix(countData = exp_data, colData = annot, design= ~status)
	    }
	  DE <- DESeq(DE)
	  name=resultsNames(DE)[2]
	  tt <- results(DE, name=name)
	}
	if(species=="mouse"){
		m2h = One2One::ortholog_data_Mouse_Human$orthologs_one2one %>% dplyr::select(human.symbol,mouse.symbol) %>% dplyr::rename(HGNC.symbol = human.symbol,MGI.symbol = mouse.symbol)
		# m2h <- readRDS("HGNC-MGI.symbols.RData")
		# Will add up this line once the code works as a package, and data can be imported(HGNC-MGI.symbols.RData)
		tt2 = cbind(tt,MGI.symbol=rownames(tt))
		tt2$MGI.symbol=as.character(tt2$MGI.symbol)
		tt3 = merge(tt2,m2h,by="MGI.symbol")
		if(!analysis_type %in% c("H_Soma-Axon","D_Soma-Axon","Soma_als-ctrl","Axon_als-ctrl")){
		  tt4 = tt3[order(tt3$P.Value),]
		} else{
		  tt4 = tt3[order(tt3$pvalue),]
		}
	}else{
		tt2 = cbind(tt,HGNC.symbol=rownames(tt))
		tt2$HGNC.symbol=as.character(tt2$HGNC.symbol)
		if(!analysis_type %in% c("H_Soma-Axon","D_Soma-Axon","Soma_als-ctrl","Axon_als-ctrl")){
		  tt4 = tt2[order(tt2$P.Value),]
		} else{
		  tt4 = tt2[order(tt2$pvalue),]
		}
	}
	write.csv(tt4,file=paste("Results_",analysis_type,"_",species,"/Tables/",file_name,sep=""))
	return(tt4)
}
