#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr desc
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom forcats fct_reorder

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


