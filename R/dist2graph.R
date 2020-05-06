
dist2graph <- function(dat.d, w=2){
  g     <- igraph::graph_from_adjacency_matrix(as.matrix(dat.d), weighted =T, mode = "upper");
  w.i   <- which(E(g)$weight > w);
  igraph::E(g)$weight[w.i] <- 0;
  return(g);
}