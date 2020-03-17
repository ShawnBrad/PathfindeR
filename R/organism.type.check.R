organism.type.check <- function(organism.type){
  if (!(organism.type %in% c('human','fly'))) stop(paste( organism.type ,' is not currently supported'))
}