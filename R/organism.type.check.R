organism.type.check <- function(organism.type){
  if (!(organism.type %in% valid.organisms)) stop(paste( organism.type ,' is not currently supported'))
}