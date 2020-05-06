require(igraph);
require(qgraph);
require(wordcloud);
## plot a nice network graph 
plot_network <- function (graph, label.gene, label.col=NULL, main=NULL){
  g2    <- graph;
  g3    <- delete_vertex_attr(g2, "name");
  e     <- get.edgelist(g3);
  
  l     <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(g3), weights=E(g2)$weight,
                                             area=8*(vcount(g3)^2),repulse.rad=(vcount(g3)^3.1))

  igraph::plot.igraph(g2, main=main, vertex.size=2, 
                      layout=l, vertex.label=NA, edge.arrow.size=0.4);
  x <- do.scale(l[,1], 1, -1);
  y <- do.scale(l[,2], 1, -1);
  cols <- rainbow(10, alpha=0.8)
  colt <- rep(cols[7], length(V(g2)$name));
  if (length(label.gene)>0){
    coli <- which(V(g2)$name %in% label.gene);
    if (length(coli)>0){
      if (!is.null(label.col)){
        colt[coli] <- label.col;
      } else {
        colt[coli] <- cols[1];
      }
    }
  }
  if (!is.na(label.gene)){
    textplot(x, y, V(g2)$name, new=F, cex=2, font=1, col=colt);
  }
}