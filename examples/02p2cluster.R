rm(list=ls());

library(gcMECM)

## load the output of Fisher's test
ct     <- "BRCA";
outf   <- paste0("clu_", ct, ".txt.gz")
inf1   <- paste0("glm_", ct, "_sum.txt.gz");
dat    <- read.table(gzfile(inf1), header=T);
## filter some genes of the genes have high frequency mutations among samples.
dat    <- do_filter(dat, g=c("TTN", "MUC16"))
## remove genes with high p values
dat.i  <- which(as.numeric(dat[,5]) < 0.001);
dat.s  <- dat[dat.i,];

## In the dat.s, column 1 and 2 were genes, and column 5 was p value by fisher's test
## convert the p values to distance matrix
dat.d  <- p2dist(dat.s[,c(1:2)], as.numeric(dat.s[,5]));
## cluster from the distance matrix 
clu    <- dist2cluster(dat.d, wt=0.00000001, method="louvain");
out.t  <- table(clu$clu$membership);
out.i  <- which(out.t>1);
out.t[out.i];
length(out.t);
## outputing the clustering results 
out.s  <- data.frame(gene=clu$clu$names, cluster=clu$clu$membership);
write.table(out.s, gzfile(outf), sep="\t", quote=F, row.names=F);
