de.vector <- function(results.table,seperate.up.down = 'none', p.cutoff = p.cutoff , lfc.cutoff = lfc.cutoff,
                      lfc.label = lfc.label , gene.label = gene.label , gene.label.type = gene.label.type,
                      padj.label = padj.label , org){

  ##### set up dbi parameters based on organism  #####
  
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

  ######## convert de genes and universe  #######
  #de.genes
  genes <- map.symbols(genes, .db = dbi ,gene.label.type = gene.label.type)

  # get gene universe and convert gene ids
  assayed.genes <- results.table %>% 
    pull(Gene) %>%
    map.symbols(.db = dbi ,gene.label.type = gene.label.type )
  
  # set up gene vector
  gene.vec=as.integer(assayed.genes %in% genes )
  names(gene.vec) = assayed.genes
  
  #gene.vec <- gene.vec[ names(gene.vec) %in% names(Org_geneIDs) ]
  
  return(gene.vec)
}

