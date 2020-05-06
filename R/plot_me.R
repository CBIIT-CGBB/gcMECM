

plot_me <- function(x, sort=T, s.lab="", main="", g.cex=1, g.col="black", m.cex=1, s.cex=1, col.bg="white", col="blue", horiz=T){
  g.col <- rep(g.col, nrow(x));
  if (sort){
    gene.i   <- order(rowSums(x), decreasing=T);
    sample.i <- rev(do.call(order, as.data.frame(t(x[gene.i,]))));
    x <- x[gene.i, sample.i];
    g.col    <- g.col[gene.i];
  }
  
  gene.n   <- nrow(x);
  sample.n <- ncol(x);
  plot(0:10, 0:10, type="n", axes=F, main=main, xlab="", ylab="");
  if (horiz){
    y.w       <- 9/gene.n;
    x.w       <- 9.5/sample.n;
    y.gap     <- y.w/20;
    gene.n1   <- gene.n+1;
    sample.n1 <- sample.n + 1;
    x.v       <- seq(1, 9.5, length.out=sample.n1);
    y.v       <- seq(1, 9, length.out=gene.n1);
    i <- 1;
    for (k in gene.n:1){
      i <- i + 1;
      j <- i - 1;
      tmp.c    <- rep(col.bg, sample.n);
      s.i      <- which(x[k,] > 0);
      tmp.c[s.i] <- col;
      xleft      <- x.v[1:sample.n];
      ybottom    <- rep(y.v[j], sample.n);
      xright     <- x.v[2:sample.n1];
      ytop       <- rep(y.v[i]-y.gap, sample.n);
      rect(xleft, ybottom, xright, ytop, col=tmp.c, border="white", lwd=0.1);
      text(0.9, y.v[j]+y.w/2, row.names(x)[k], pos=2, cex=g.cex, col=g.col[k]);
    }
    text(mean(x.v), 0.5, s.lab, cex=s.cex);
    text(mean(x.v), 9.5, main, cex=m.cex);
  } else {
    x.w       <- 9/gene.n;
    y.w       <- 9.5/sample.n;
    x.gap     <- x.w/20;
    gene.n1   <- gene.n+1;
    sample.n1 <- sample.n + 1;
    x.v       <- seq(1, 9, length.out=gene.n1);
    y.v       <- rev(seq(1, 9, length.out=sample.n1));
    i <- 1;
    for (k in 1:gene.n){
      i <- i + 1;
      j <- i - 1;
      tmp.c    <- rep(col.bg, sample.n);
      s.i      <- which(x[k,] > 0);
      tmp.c[s.i] <- col;
      xleft      <- rep(x.v[j], sample.n);
      ybottom    <- rep(y.v[1:sample.n]);
      xright     <- rep(x.v[i], sample.n);
      ytop       <- rep(y.v[2:sample.n1]);
      rect(xleft, ybottom, xright, ytop, col=tmp.c, border="white", lwd=0.1);
      text(x.v[i]-x.w/4, 0.9, row.names(x)[k], pos=2, cex=g.cex, srt=90, col=g.col[k]);
    }
    text(0.5, mean(y.v), s.lab, cex=s.cex, srt=90);
    text(mean(x.v), 9.5, main, cex=m.cex);
  }
}