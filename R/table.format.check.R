table.format.check <- function(results.table, gene.label){
  if (!(is.data.frame(results.table))) stop('Entered results table is not a data frame', call. = F)
  if (!( gene.label %in% colnames(results.table))) stop(paste(gene.label,
                                                              " column not present in table, check gene.label argument" ,call. = F))
}

