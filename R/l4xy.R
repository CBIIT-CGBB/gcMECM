

do.scale <- function(x, newMax, newMin){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin};

## l: layouts; xy: coordination; r:radius
## output: the layouts will be in the region: centre=xy and radius=r
l4xy <- function(l, xy, r){
  newMax  <- xy[1]+r[1];
  newMin  <- xy[1]-r[1];
  x       <- do.scale(l[,1], newMax, newMin);
  newMax  <- xy[2]+r[2];
  newMin  <- xy[2]-r[2];
  y       <- do.scale(l[,2], newMax, newMin);
  out     <- data.frame(x=x, y=y);
  out     <- as.matrix(out);
  return(out);
}