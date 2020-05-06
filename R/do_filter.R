
do_filter <- function(dat, g){
  i <- which(dat[,1] %in% g);
  j <- which(dat[,2] %in% g);
  k <- unique(c(i, j));
  if (length(k)>0){
    dat <- dat[-k,];
  } else {
    dat <- dat;
  }
  return(dat);
}