#' Plot Summary of GSA adjusting the threshold P-value using Bonferroni correction.
#'
#' \code{plot_GSA} Creates a barplot of -log(P-value) of the GSA analysis for 250,300 and 400 DR and UR gene sets.
#'
#' @param summary_filepath file path including file name of the summary file containing at all the gsa.out results. At least the P-value, NGENES, and variable name gsa.out must be present. The function assumes that the file does not contain comments.
#' @param sep separator in the file, can be set to tab="\t", space=" ", etc. Default is tab.
#' @param plot_name the name of the jpeg/pdf files to be created, without the extension.
#' @param plot_maintitle the main title of the plots.
#' @param width width of the jpeg file. Default is 1000.
#' @param height height of the jpeg file. Default is 800.
#'
#' @return One pdf and one jpeg plot files.
#'
#' @examples
#' # Plots the tissue analysis results considering a default tab separator and default height and width.
#' plot_GSA(file_path="/Users/bioinfo/all_results.txt", plot_name="GSA_Mouse_HSomaAxon",plot_maintitle="GSA Mouse Healthy Soma-Axon")
#'
plot_GSA <- function(summary_filepath,sep='\t',plot_name='GSA_Results',plot_maintitle='GSA Results',width=1000,height=800) {
  GSA <- read.csv(summary_filepath,sep=sep,row.names=1)
  #Read file and convert it into a data frame whose row.names are equal to the gene symbols.
  logneg <- -log(GSA$P)
  GSA['-logP']=logneg
  #Calculates the -log of the GSA P-values and creates a new variable in the data frame.
  file_location <- dirname(dirname(summary_filepath))
  if (!file.exists(paste(file_location,"/Figures/",sep = ""))){
    dir.create(paste(file_location,"/Figures/",sep = ""))
  }
  jpeg(paste(file_location,'/Figures/',plot_name,".jpeg",sep=""), width = width, height = height)
  op <- par(mar=c(9,6,4,4))
  col_list=c(replicate(3,'red'),replicate(3,'steelblue'))
  y_lim <- -log(0.05/400)+2
  if(max(GSA$`-logP`) > y_lim){
    y_lim = max(GSA$`-logP`)+2
  }
  bar=barplot(GSA$`-logP`,names.arg = row.names(GSA),las=2,ylab = '-log(Pval)',ylim=c(0,y_lim),col=col_list,main=plot_maintitle)
  #legend('right', legend=c('DR','UR'), fill=c('red','steelblue'),cex=0.75)
  text(bar, GSA$`-logP`, labels=GSA$NGENES, pos=3)
  # Adds the number of genes in the set to the bar.
  legend(title="Bonf. sign.",'right', legend=c('Nominal','250','300','400'), col=c('black','red','blue','orange'),lty=2,cex=1)
  abline(h=-log(0.05),col='black',lty='dotted',lwd=3)
  abline(h=-log(0.05/250),col='red',lty='dotted',lwd=3)
  abline(h=-log(0.05/300),col='blue',lty='dotted',lwd=3)
  abline(h=-log(0.05/400),col='orange',lty='dotted',lwd=3)
  dev.off()
  pdf(paste(file_location,"/Figures/",plot_name,".pdf",sep=""),width=11)
  op <- par(mar=c(8,4,4,4))
  col_list=c(replicate(3,'red'),replicate(3,'steelblue'))
  bar=barplot(GSA$`-logP`,names.arg = row.names(GSA),las=2,ylab = '-log(Pval)',ylim=c(0,y_lim),col=col_list,main=plot_maintitle)
  #legend('right', legend=c('DR','UR'), fill=c('red','steelblue'),cex=0.75)
  text(bar, GSA$`-logP`, labels=GSA$NGENES, pos=3)
  # Adds the number of genes in the set to the bar.
  legend(title="Bonf. sign.",'right', legend=c('Nominal','250','300','400'), col=c('black','red','blue','orange'),lty=2,cex=1)
  abline(h=-log(0.05),col='black',lty='dotted',lwd=3)
  abline(h=-log(0.05/250),col='red',lty='dotted',lwd=3)
  abline(h=-log(0.05/300),col='blue',lty='dotted',lwd=3)
  abline(h=-log(0.05/400),col='orange',lty='dotted',lwd=3)
  dev.off()
}



