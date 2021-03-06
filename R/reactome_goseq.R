#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import goseq
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom stats p.adjust
#' @importFrom AnnotationDbi get

reactome_goseq <- function(gene.vector, pwf, org = organism.type, use.method){


  if (org == 'fly'){

    Pathways <- DM_pathways
    Org_geneIDs = DM_geneIDs
    pathName_prefix <- 'Drosophila melanogaster: '
    org.db <- org.Dm.egSYMBOL
    dbi = org.Dm.eg.db
  }


  if (org == 'human'){
   
    Pathways <- HS_pathways
    Org_geneIDs = HS_geneIDs
    pathName_prefix <- 'Homo sapiens: '
    org.db <- org.Hs.egSYMBOL
    dbi = org.Hs.eg.db
  }
  
  # run go seq
  temp.res <- suppressMessages(goseq(pwf, gene2cat = Org_geneIDs, method = use.method)) %>%
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
    allegs<- suppressMessages(AnnotationDbi::get(.term, reactomePATHID2EXTID)) %>%
      entrez.to.symbols(dbi)
    
    term.cg <- intersect(allegs,de.genes)
    out.genes <- c(out.genes, toString(term.cg))
  }
  temp.res$Genes <- out.genes
  return(temp.res)
}







