pwf <- function(gene.vector, org = organism.type){

  if (org == 'fly') ( genome = 'dm3')
  if (org == 'human') ( genome = 'hg19')

  pwf = nullp(gene.vector,genome = genome, id = 'geneSymbol' ,plot.fit = F) %>%
    return()
}

