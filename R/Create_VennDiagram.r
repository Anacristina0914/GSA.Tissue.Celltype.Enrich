#' Creates and stores in a current directory a Venn Diagram using 2 gene/element lists.
#'
#' \code{Create_VennDiagram} Creates a Venn Diagram using two lists of elements (genes,attributes,cats,pokemons,etc).
#'
#' @param GeneList1 A list of elements.
#' @param GeneList2 A list of elements.
#' @param varname_1 A string containing the variable name of the first variable to be displayed in the plot.
#' @param varname_2 A string containing the variable name of the second variable to be displayed in the plot.
#' 
#' @return A 350x350 png image of a Venn Diagram of the elements in list1 and list2. Image is saved in the current working directory.
#'
#' @examples
#' list1 <- c("dog","cat","mouse","ostrich","mole")
#' list2 <- c("mole","dog","bird","shark")
#' Create_VennDiagram(list1,list2,"Animals1","Animals2")
#'

Create_VennDiagram <- function(GeneList1,GeneList2,varname_1,varname_2){
  if (!require("VennDiagram")){
    install.packages("VennDiagram")
  }
  if (!require("RColorBrewer")){
    install.packages("RColorBrewer")
  }
  library(VennDiagram)
  library(RColorBrewer)
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  #Avoids the creation of a log file everytime a new VennDiagram is created.
  venn.diagram(
    x=list(GeneList1,GeneList2),
    category.names = c(varname_1,varname_2),
    filename=paste("Venn_",varname_1,varname_2,".png",sep=""),
    #Graphic parameters
    height = 350,
    width= 350,
    lwd = 2,
    lty = 'blank',
    fill = c("#B3E2CD","#FDCDAC"),
    #Numbers
    cex = 0.18,
    #Names
    cat.cex = 0.18,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    cat.pos=c(-27,27)
  )
}
