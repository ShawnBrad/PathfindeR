#' @export
reactome_lolliplot <- function(res.reactome){
  colp <- c(up = 'darkred', down = 'darkblue')
  res.reactome %>%
    mutate(Ratio = numDEInCat/numInCat) %>%
    mutate(padj = log10(padj)) %>%
    mutate(scaling = if_else(FC.direction == 'down',1,-1)) %>%
    mutate(padj = padj * scaling) %>%
    arrange(desc(padj)) %>%
    
    ggplot(aes(fct_reorder(path_name, padj), padj, size = Ratio, color = FC.direction)) +
    geom_segment(aes(x=fct_reorder(path_name, padj), xend=fct_reorder(path_name, padj),
                     y=0, yend=padj), color="grey25", size = .5) +
    ylab('log10 FDR')+ xlab('Enriched Pathway') +
    scale_color_manual(name = 'Up vs Down Regulated',values = colp)+
    geom_point( ) +
    coord_flip() %>%
    return()
}