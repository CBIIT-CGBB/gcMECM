rm(list=ls());

rm(list=ls());

library(gcMECM);
options(stringsAsFactors = F);

f.s  <- "BRCA";
inf1 <- paste0("clu_", f.s, ".txt.gz");
dat1 <- read.table(inf1, header=T);
dat1.t <- table(dat1[,2]);
dat1.i <- which(dat1.t > 50);
clu1.i <- names(dat1.t)[dat1.i];

inf2   <- paste0("mutect_", f.s, "_matrix.txt.gz");
dat2   <- read.table(inf2, header=T, row.names=1);

pdff   <- paste0("04_1plot_ME.pdf");
pdf(pdff, 18, 12);
for (cl.i in clu1.i[3]){
  i      <- which(dat1[,2]==cl.i);
  gene   <- dat1[i,1];
  dat2.i <- which(row.names(dat2) %in% gene);
  dat2.s <- dat2[dat2.i,];
  dat2.s <- dat2.s[1:80,]
  dat2.n <- apply(dat2.s, 2, sum);
  dat2.j <- which(dat2.n==0);
  dat2.s <- dat2.s[,-dat2.j];
  gene.s <- length(dat2.i);
  samp.s <- ncol(dat2.s);
  main   <- paste0(f.s, " clusterID:", cl.i, " Gene#:", gene.s, " sample#:", samp.s);
  plot_me(dat2.s, g.cex=0.2, g.col="black", col.bg="gray90", col="blue", main=main);
}
dev.off();

