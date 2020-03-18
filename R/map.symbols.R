#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom magrittr %>%

map.symbols<- function(genes, .db, gene.label.type){
  ##### set up dbi parameters based on organism  #####
  
  
  genes <- suppressMessages(mapIds(.db, keys = genes,
                                 column = "SYMBOL", keytype = gene.label.type, multiVals = 'first'))  %>%
    as.character()
  # remove nas
  genes<- genes[!is.na(genes)]
  genes<- genes[!duplicated(genes)]
  return(genes)
  
}

