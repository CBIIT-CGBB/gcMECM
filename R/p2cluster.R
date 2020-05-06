require(igraph);

p2cluster <- function(dat=dat, method="louvain", wt=NULL){
  g     <- graph_from_data_frame(dat, directed = F, vertices = NULL);
  w.j   <- which(is.na(E(g)$weight));
  if (length(w.j) > 0){
    g       <- delete.edges(g, w.j);
  }
  if (!is.null(wt)){
    w.i   <- which(E(g)$weight > wt);
    if (length(wi)>0){
      igraph::E(g)$weight[w.i] <- 0;
    }
  }
  if (method=="louvain"){
    out     <- igraph::cluster_louvain(g, weights = E(g)$weight);  
  } else if (method=="fast_greedy"){
    out     <- igraph::cluster_fast_greedy(g, weights = E(g)$weight);
  } else if (method=="infomap"){
    out     <- igraph::cluster_infomap(g, e.weights = E(g)$weight);
  } else if (method=="label_prop"){
    out     <- igraph::cluster_label_prop(g, weights = E(g)$weight); 
  } else if (method=="spinglass"){
    out     <- igraph::cluster_spinglass(g, weights = E(g)$weight); 
  } else if (method=="weight"){
    g       <- delete.edges(g, which(E(g)$weight == 0))
    out     <- igraph::components(g); 
  } else {
    stop("The method is louvain, fast_greedy, infomap, label_prop or spinglass.\n")
  }
  return(list(clu=out, g=g));
}