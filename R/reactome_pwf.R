pwf <- function(gene.vector, id = 'geneSymbol', org){

  if (org == 'fly') ( genome = 'dm3')
  if (org == 'human') ( genome = 'hg19')

  pwf = nullp(gene.vector,genome = genome, id = id ,plot.fit = F) %>%
    return()
}

