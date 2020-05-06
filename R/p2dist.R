p2dist <- function(dat=dat, p=p){
  ## check if P value = 0;
  p.i    <- which(p==0);
  if (length(p.i) > 0){
    p.min <- min(p[-p.i], na.rm =TRUE);
    p[p.i] <- p.min;
  }
  w     <- p;
  dat.s <- data.frame(dat[,c(1,2)], w=w);
  dat.s <- rbind(dat.s, dat.s[,c(2,1,3)]);
  n.u   <- unique(c(as.character(dat.s[,1]), as.character(dat.s[,2])));
  n.n   <- length(n.u);
  dat.m <- matrix(rep(NA, n.n*n.n), ncol=n.n);
  colnames(dat.m) <- n.u;
  rownames(dat.m) <- n.u;
  
  col.i <- match(dat.s[,1], n.u);
  row.i <- match(dat.s[,2], n.u);
  
  for (i in 1:length(col.i)){
    i1 <- col.i[i]
    i2 <- row.i[i];
    dat.m[i1, i2] <- dat.s[i,3];
  }
  diag(dat.m) <- 0;
  dat.m[lower.tri(dat.m)] <- t(dat.m)[lower.tri(dat.m)];
  dat.d <- as.dist(dat.m);
  return(dat.d);
}