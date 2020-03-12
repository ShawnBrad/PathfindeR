de.vector <- function(results.table,seperate.up.down = 'none', p.cutoff , lfc.cutoff  ,
                      lfc.label , gene.label , gene.label.type ,
                      padj.label , org){


if (org == 'fly') (dbi = drosophila2.db)
if (org == 'human') (dbi = org.Hs.eg.db)


  lfc.label <- as.symbol(lfc.label)
  gene.label <- as.symbol(gene.label)
  padj.label <- as.symbol(padj.label)

  # get de genes from results

  if (seperate.up.down == 'up') {
    genes <- results.table %>%
      filter(!!padj.label <= p.cutoff) %>%
      filter(!!lfc.label >= lfc.cutoff) %>%
      pull(!!gene.label)

  } else if (seperate.up.down == 'down') {
      genes <- results.table %>%
        filter(!!padj.label <= p.cutoff) %>%
        filter(!!lfc.label <= -1*lfc.cutoff) %>%
        pull(!!gene.label)

  } else if (seperate.up.down == 'none'){
    genes <- results.table %>%
      filter(!!padj.label <= p.cutoff) %>%
      filter(abs(!!lfc.label) >= lfc.cutoff) %>%
      pull(!!gene.label)
  }



  genes <- AnnotationDbi::select(dbi, keys = genes,
                                 columns = "SYMBOL", keytype = gene.label.type)  %>%
    dplyr::pull(SYMBOL)

  # remove nas
  genes<- genes[!is.na(genes)]

  # get gene universe
  assayed.genes <- results.table %>%
    pull(Gene)

  assayed.genes <- AnnotationDbi::select(dbi, keys = assayed.genes,
                                         columns = "SYMBOL", keytype = gene.label.type)  %>%
    dplyr::pull(SYMBOL)

  # remove nas
  assayed.genes<- assayed.genes[! is.na(assayed.genes)]

  # set up gene vector
  gene.vec=as.integer(assayed.genes %in% genes )
  names(gene.vec) = assayed.genes
  return(gene.vec)
}

