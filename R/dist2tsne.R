
dist2tsne <- function(dat.d, tsne.per, tsne.max){
  tsne  <- Rtsne::Rtsne(dat.d, is_distance =T,  perplexity = tsne.per, max_iter = tsne.max);
  return(tsne);
}
