#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr desc
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 geom_segment
#' @importFrom forcats fct_reorder
#' @importFrom base return
#' @importFrom base log10

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