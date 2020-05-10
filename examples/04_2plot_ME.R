rm(list=ls());

library(gcMECM);
options(stringsAsFactors = F);

f.s  <- "BRCA";
pdff <- paste0("04_2plot_ME.pdf");
pdf(pdff, 5, 12);
par(mfrow=c(3, 1), mar=c(1,2,1,1));
## the file was from 03cluster_map
infile <- paste0("BRCA.txt")
dat    <- read.table(infile, header=T);
clu.s  <- paste0(dat[,1], ".", dat[,2]);
clu.t  <- table(clu.s);
clu.i  <- which(clu.t > 1);
clu.u  <- names(clu.t[clu.i]);
mu.f   <- paste0("mutect_", f.s, "_matrix.txt.gz");
mdat   <- read.table(gzfile(mu.f), header=T, row.names=1);
mdat[mdat>1] <- 1;
for (clu.x in clu.u[c(2,3,6)]){
  clu.j  <- which(clu.s==clu.x);
  f.gene <- dat[clu.j, 3];
  mdat.i <- which(row.names(mdat) %in% f.gene);
  mdat.s <- mdat[mdat.i,];
  coln   <- colSums(mdat.s);
  coln.i <- which(coln==0);
  mdat.s <- mdat.s[,-coln.i];
  main   <- paste(f.s, clu.x);
  plot_me(mdat.s, g.cex=0.6, g.col="black", col.bg="gray90", col="blue", main="");
}
dev.off();
