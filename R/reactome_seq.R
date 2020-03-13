

reactome_seq <- function(results.table, seperate.up.down = TRUE, p.cutoff = 0.05, lfc.cutoff = 0.7 ,
                         lfc.label = 'log2FoldChange', gene.label = 'Gene', gene.label.type = "FLYBASE",
                         padj.label = 'padj', organism.type = 'fly'){

  lfc.label <- as.symbol(lfc.label)
  gene.label <- as.symbol(gene.label)
  padj.label <- as.symbol(padj.label)


  if (seperate.up.down){
    print('separating by up vs down fold change')

    gene.vector <- de.vector(results.table,seperate.up.down = 'up', p.cutoff , lfc.cutoff ,
                             lfc.label, gene.label, gene.label.type, padj.label, org = organism.type)

    pwf.model <- pwf (gene.vector,org = organism.type)
    reactome.results.up <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type)
    #print(reactome.results.up)


    # reactome for down genes
    gene.vector <- de.vector(results.table,seperate.up.down = 'down', p.cutoff , lfc.cutoff ,
                             lfc.label, gene.label, gene.label.type, padj.label, org = organism.type)

    pwf.model <- pwf (gene.vector,org = organism.type)
    reactome.results.down <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type )

    #print(reactome.results.down)

    reactome.results <- reactome.results.up %>%
      add_column(FC.direction = 'up') %>%
      bind_rows( reactome.results.down %>%
                   add_column(FC.direction = 'down'))


  }
  else {

    gene.vector <- de.vector(results.in,seperate.up.down = 'none',  org = organism.type)

    pwf.model <- pwf (gene.vector,org = org)
    reactome.results <- reactome_goseq(gene.vector, pwf = pwf.model, org = organism.type)

  }

  return(reactome.results)
  }
