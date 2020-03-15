entrez.to.symbols<- function(genes, .db){
  mapIds(.db, keys = genes,
                        column = "SYMBOL", keytype = "ENTREZID", multiVals = 'first')  %>%
    as.character() %>%
    return()
}

