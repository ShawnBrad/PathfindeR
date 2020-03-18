#' @export
reactome_barplot <- function(res.reactome){
  colp <- c(up = 'darkred', down = 'darkblue')
  res.reactome %>%
    mutate(padj = log10(padj)) %>%
    mutate(scaling = if_else(FC.direction == 'down',1,-1)) %>%
    mutate(padj = padj * scaling) %>%
    arrange(desc(padj)) %>%
    
    ggplot(aes(fct_reorder(path_name, padj), padj, fill = FC.direction)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(name = 'Up vs Down Regulated',values = colp)+
    ylab('log10 FDR')+ xlab('Enriched Pathway') +
    coord_flip() %>%
    return()
}

