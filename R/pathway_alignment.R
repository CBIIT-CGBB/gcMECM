
require(igraph);
pathway_alignment <- function(g0, relations){
  g      <- graph_from_data_frame(relations, directed = TRUE, vertices = NULL);
  g.j    <- which(V(g)$name %in% genes);
  g1     <- (g0 %s% g)
  g1.c   <- components(g1);
  g1.t   <- table(g1.c$membership);
  g1.i   <- which(g1.t > 2);
  if (length(g1.i) < 1){
    return(NA);
  }
  ## mapped genes
  out.n <- NULL;
  ## relationships of mapped genes
  out.s <- NULL;
  for (c.i in names(g1.i)){
    clu.j <- which(g1.c$membership==as.numeric(c.i));
    g2    <- induced.subgraph(g1, clu.j);
    g.n <- length(V(g2)$name);
    df.n <- as_data_frame(g2);
    df.m <- data.frame(sub=rep(c.i, nrow(df.n)),
                       df.n);
    df.s <- data.frame(sub=rep(c.i, g.n),
                       gene=V(g2)$name);
    out.s <- rbind(out.s, df.s);
    out.n <- rbind(out.n, df.m);
  }
  return(list(out.s, out.n))
}