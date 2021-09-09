#' Plot Tissue Enrichment Analysis results adjusting P-value using FDR.
#'
#' \code{plot_MAGMATissueAnalysis_FDR} Saves a barplot using P-value by tissue corrected using FDR.
#'
#' @param file_path file path includying file name of the gsa.out tissue enrichment file. The function assumes that the file does not contain commnents, is tab-separated and tissues are ordered in alphabetical order, considering the brain areas first and the non-brain areas last.
#' @param sep separator in the file, can be set to tab="\t", space=" ", etc.
#' @param plot_name the name of the jpeg/pdf files to be created, without the extension.
#' @param plot_maintitle the main title of the plots.
#' @param width width of the jpeg file.
#' @param height height of the jpeg file.
#'
#' @return two plot files, one jpeg and one pdf files.
#'
#' @examples
#' # Plots the tissue analysis results considering a default tab separator and default height and width.
#' plot_MAGMATissueAnalysis_FDR(file_path="/Users/bioinfo/tissue_analysis.gsa.out",plot_name="TissueEnrich_Results_FDR",plot_maintitle="Tissue Enrichment Analysis MAGMA FDR Corrected")
#' 
plot_MAGMATissueAnalysis_FDR <- function(file_path,sep='\t',plot_name='TissueEnrich_Results_FDR.jpeg',plot_title='Tissue Enrichment Analysis MAGMA FDR Corrected MAF > X',width=1000,height=800) {
	Tissue_enrich <- read.csv(file_path,sep='\t',row.names=1)
	#Read file and convert it into a dataframe whose row.names are equal to the gene symbols.
	pval_nom_BH <- -log(p.adjust(Tissue_enrich$P,method='fdr',n=length(Tissue_enrich$P)))
	#Calculates -log of FDR from P-val.
	Tissue_enrich['-log(FDR)']=pval_nom_BH
	file_location <- dirname(file_path)
	if (!file.exists(paste(file_location,"/Graphs/",sep = ""))){
		dir.create(paste(file_location,"/Graphs/",sep = ""))
	}
	jpeg(paste(file_location,'/Graphs/',plot_name,sep=""), width = width, height = height)
	op <- par(mar=c(18,6,4,4))
	col_list=c(replicate(13,'steelblue'),replicate(24,'coral'))
	barplot(Tissue_enrich$`-log(FDR)`,names.arg = Tissue_enrich$FULL_NAME,las=2,ylab = '-log(FDR)',ylim=c(0,5),col=col_list,main=plot_title)
	legend('right', legend=c('CNS','Non-CNS'), fill=c('steelblue','coral'),title='Tissue Localization',cex=0.75)
	abline(h=-log(0.05),col='red',lty='dotted')
	dev.off()
	pdf(paste(file_location,"/Graphs/",plot_name,".pdf",sep=""),width=11)
	op <- par(mar=c(17,4,4,4))
	col_list=c(replicate(13,'steelblue'),replicate(24,'coral'))
	barplot(Tissue_enrich$`-logP`,names.arg = Tissue_enrich$FULL_NAME,las=2,ylab = '-log(Pval)',ylim=c(0,8),col=col_list,main=plot_title)
	legend('right', legend=c('CNS','Non-CNS'), fill=c('steelblue','coral'),title='Tissue Localization',cex=0.75)
	abline(h=pval_nom_bonf,col='red',lty='dotted')
	dev.off()
}
