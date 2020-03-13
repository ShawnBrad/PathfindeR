de.vector <- function(results.table,seperate.up.down = 'none', p.cutoff = p.cutoff , lfc.cutoff = lfc.cutoff  ,
                      lfc.label = lfc.label , gene.label = gene.label , gene.label.type =gene.label.type ,
                      padj.label = padj.label , org = organism.type){

  
  if (org == 'fly'){
    load("data/DM_pathways.RData")
    load("data/DM_geneIDs.RData")
    Pathways <- DM_pathways
    Org_geneIDs = DM_geneIDs
    pathName_prefix <- 'Drosophila melanogaster: '
    org.db <- org.Dm.egSYMBOL
    dbi = drosophila2.db
  }
  
  
  if (org == 'human'){
    load("data/HS_pathways.RData")
    load("data/HS_geneIDs.RData")
    Pathways <- HS_pathways
    Org_geneIDs = HS_geneIDs
    pathName_prefix <- 'Homo sapiens: '
    org.db <- org.Hs.egSYMBOL
    dbi = org.Hs.eg.db
  }
  



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



  genes <- AnnotationDbi::mapIds(dbi, keys = genes,
                                 column = "SYMBOL", keytype = gene.label.type, multiVals = 'first')  %>%
    as.character()

  # remove nas
  genes<- genes[!is.na(genes)]
  genes<- genes[!duplicated(genes)]
  
  # get gene universe
  assayed.genes <- results.table %>%
    pull(Gene)

  assayed.genes <- AnnotationDbi::mapIds(dbi, keys = assayed.genes,
                                         column = "SYMBOL", keytype = gene.label.type, multiVals = 'first')  %>%
    as.character()

  # remove nas
  assayed.genes<- assayed.genes[! is.na(assayed.genes)]
  assayed.genes<- assayed.genes[! duplicated(assayed.genes)]
  # set up gene vector
  gene.vec=as.integer(assayed.genes %in% genes )
  names(gene.vec) = assayed.genes
  
  gene.vec <- gene.vec[ names(gene.vec) %in% names(Org_geneIDs) ]
  
  return(gene.vec)
}

