#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom magrittr %>%

entrez.to.symbols<- function(genes, .db){
  suppressMessages(mapIds(.db, keys = genes,
                        column = "SYMBOL", keytype = "ENTREZID", multiVals = 'first'))  %>%
    as.character() %>%
    return()
}
