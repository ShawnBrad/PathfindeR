reactome_goseq <- function(gene.vector, pwf, org = 'fly' ){


  if (org == 'fly'){
    load("data/DM_pathways.RData")
    load("data/DM_geneIDs.RData")
    Pathways <- DM_pathways
    Org_geneIDs = DM_geneIDs
    pathName_prefix <- 'Drosophila melanogaster: '
    org.db <- org.Dm.egSYMBOL
  }


  if (org == 'human'){
    load("data/HS_pathways.RData")
    load("data/HS_geneIDs.RData")
    Pathways <- HS_pathways
    Org_geneIDs = HS_geneIDs
    pathName_prefix <- 'Homo sapiens: '
    org.db <- org.Hs.egSYMBOL
  }


  # run go seq
  temp.res <- goseq(pwf, gene2cat = Org_geneIDs) %>%
    mutate(padj= p.adjust(over_represented_pvalue, method="BH"))  %>%
    filter(padj <= 0.05) %>%
    left_join(Pathways,  by = c('category' = 'DB_ID')) %>%
    mutate(path_name = str_remove(path_name,pathName_prefix ))%>%
    dplyr::select( - over_represented_pvalue, -under_represented_pvalue )


  ## get list of de genes
  de.genes <- gene.vector[gene.vector == 1] %>%
    names()

  ## get genes from cats
  out.genes<- c()
  .id<- temp.res$category

  # get gene ist for each category
  for (.term in .id){
    allegs<-AnnotationDbi::get(.term, reactomePATHID2EXTID)
    genes = unlist(mget(allegs,org.db),use.names = FALSE)
    term.cg <- intersect(genes,de.genes)
    out.genes <- c(out.genes, toString(term.cg))

  }
  temp.res$Genes <- out.genes
  return(temp.res)
}







