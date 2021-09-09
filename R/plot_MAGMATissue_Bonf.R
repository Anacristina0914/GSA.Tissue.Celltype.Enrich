#' Plot Tissue Enrichment Analysis results adjusting P-value using Bonferroni correction.
#'
#' \code{plot_MAGMATissueAnalysis_Bonf} Saves a barplot using P-value by tissue corrected using Bonferroni correction.
#'
#' @param file_path file path includying file name of the gsa.out tissue enrichment file. The function assumes that the file does not contain commnents and is tab-separated.
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
#' plot_MAGMATissueAnalysis_Bonf(file_path="/Users/bioinfo/tissue_analysis.gsa.out",plot_name="TissueEnrich_Results_Bonf",plot_maintitle="Tissue Enrichment Analysis MAGMA Bonf. Corrected")
plot_MAGMATissueAnalysis_Bonf <- function(file_path,sep='\t',plot_name='TissueEnrich_Results_Bonf',plot_maintitle='Tissue Enrichment Analysis MAGMA Bonf. Corrected',width=1000,height=800) {
	Tissue_enrich <- read.csv(file_path,sep='\t',row.names=1)
	#Read file and convert it into a dataframe whose row.names are equal to the gene symbols.
	logneg <- -log(Tissue_enrich$P)
	Tissue_enrich['-logP']=logneg
	#Calculates the -log of the Tissue Prioritization Pvalues and creates a new variable in the data frame.
	pval_nom_bonf <- -log(0.05/length(row.names(Tissue_enrich)))
	#Calculates P-val significance in log scale correcting for the number of Tissues in the analysis.
	file_location <- dirname(file_path)
	if (!file.exists(paste(file_location,"/Graphs/",sep = ""))){
		dir.create(paste(file_location,"/Graphs/",sep = ""))
	}
	jpeg(paste(file_location,'/Graphs/',plot_name,".jpeg",sep=""), width = width, height = height)
	op <- par(mar=c(18,6,4,4))
	col_list=c(replicate(13,'steelblue'),replicate(24,'coral'))
	barplot(Tissue_enrich$`-logP`,names.arg = Tissue_enrich$FULL_NAME,las=2,ylab = '-log(Pval)',ylim=c(0,8),col=col_list,main=plot_title)
	legend('right', legend=c('CNS','Non-CNS'), fill=c('steelblue','coral'),title='Tissue Localization',cex=0.75)
	abline(h=pval_nom_bonf,col='red',lty='dotted')
	dev.off()
	pdf(paste(file_location,"/Graphs/",plot_name,".pdf",sep=""),width=11)
	op <- par(mar=c(17,4,4,4))
	col_list=c(replicate(13,'steelblue'),replicate(24,'coral'))
	barplot(Tissue_enrich$`-logP`,names.arg = Tissue_enrich$FULL_NAME,las=2,ylab = '-log(Pval)',ylim=c(0,8),col=col_list,main=plot_title)
	legend('right', legend=c('CNS','Non-CNS'), fill=c('steelblue','coral'),title='Tissue Localization',cex=0.75)
	abline(h=pval_nom_bonf,col='red',lty='dotted')
	dev.off()
}
