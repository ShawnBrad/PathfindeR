#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import goseq
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom dplyr add_columns




#' @export
reactome_seq <- function(results.table, seperate.up.down = TRUE, p.cutoff = 0.05, lfc.cutoff = 0.7 ,
                         lfc.label = 'log2FoldChange', gene.label = 'Gene', gene.label.type = "ENSEMBL",
                         padj.label = 'padj', organism.type = 'human', goseq.method = "Wallenius"){
  
  ## Check for valid parameters 
  organism.type.check(organism.type) 
  table.format.check(results.table , gene.label)
  gene.label.check(gene.label.type, organism.type)
  
  
  ##### set up dbi parameters based on organism  #####
  
  if (organism.type == 'fly'){
    
    Pathways <- DM_pathways
    Org_geneIDs = DM_geneIDs
    pathName_prefix <- 'Drosophila melanogaster: '
    org.db <- org.Dm.egSYMBOL
    dbi = org.Dm.eg.db
  }
  
  
  if (organism.type == 'human'){
    
    Pathways <- HS_pathways
    Org_geneIDs = HS_geneIDs
    pathName_prefix <- 'Homo sapiens: '
    org.db <- org.Hs.egSYMBOL
    dbi = org.Hs.eg.db
  }
  #######

  lfc.label <- as.symbol(lfc.label)
  gene.label <- as.symbol(gene.label)
  padj.label <- as.symbol(padj.label)


  if (seperate.up.down){
    print('separating by up vs down fold change')

    gene.vector <- de.vector(results.table,seperate.up.down = 'up',p.cutoff =  p.cutoff , lfc.cutoff = lfc.cutoff ,
                                  lfc.label = lfc.label, gene.label = gene.label, gene.label.type = gene.label.type, 
                             padj.label = padj.label, org = organism.type)

    pwf.model <- pwf (gene.vector,org = organism.type)
    reactome.results.up <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type, use.method = goseq.method)
   # print('finished up results' )


    # reactome for down genes
    gene.vector <- de.vector(results.table,seperate.up.down = 'down',p.cutoff =  p.cutoff , lfc.cutoff = lfc.cutoff ,
                             lfc.label = lfc.label, gene.label = gene.label, gene.label.type = gene.label.type, 
                             padj.label = padj.label, org = organism.type)

    pwf.model <- pwf (gene.vector,org = organism.type)
    reactome.results.down <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type ,use.method = goseq.method)

    #print('finished down results' )

    reactome.results <- reactome.results.up %>%
      add_column(FC.direction = 'up') %>%
      bind_rows( reactome.results.down %>%
                   add_column(FC.direction = 'down'))


  }
  else {

    gene.vector <- de.vector(results.table,seperate.up.down = 'none',p.cutoff =  p.cutoff , lfc.cutoff = lfc.cutoff ,
                             lfc.label = lfc.label, gene.label = gene.label, gene.label.type = gene.label.type, padj.label = padj.label, org = organism.type)

    pwf.model <- pwf (gene.vector,org =  organism.type)
    reactome.results <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type, use.method = goseq.method)

  }

  return(reactome.results)
  }
