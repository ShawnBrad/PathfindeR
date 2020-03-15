gene.label.check <- function(gene.label.type, organism.type){
  # Check for support if fly
  if (organism.type == 'fly'){
    if (!(gene.label.type %in% keytypes(org.Dm.eg.db))) stop (paste(gene.label.type,
      'is not a valid gene label type, use keytypes(org.Dm.eg.db) to find a list of valid identifiers'),call. = F)
  }
  # Check for support if human
  if (organism.type == 'human'){
    if (!(gene.label.type %in% keytypes(org.Hs.eg.db))) stop (paste(gene.label.type,
     'is not a valid gene label type, use keytypes(org.Hs.eg.db) to find a list of valid identifiers'),call. = F)
  }
}




