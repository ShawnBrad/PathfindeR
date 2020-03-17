reactome_bubbleplot <- function(res.reactome){
  
  res.reactome %>%
    mutate(Ratio = numDEInCat/numInCat) %>%
    mutate(padj = log10(padj)) %>%
    mutate(scaling = if_else(FC.direction == 'down',1,-1)) %>%
    mutate(padj = padj * scaling) %>%
    arrange(desc(padj)) %>%
    
    ggplot(aes(Ratio, padj, size = numDEInCat, color = path_name, label = path_name)) +
    geom_point(alpha = .6,show.legend = F)+
    scale_size(range = c(4,10)) +
    geom_text(check_overlap = T,nudge_y = 1.2, show.legend = F, size = 3, color = 'grey25')+
    ylab('log10 FDR')+ xlab('Ration of DE genes to genes in category') %>%
    return()
}



