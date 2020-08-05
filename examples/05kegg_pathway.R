rm(list=ls());

options(stringsAsFactors = F);

library(gcMECM);
library(OmicPath);
library(NCIRASPathway);

g.dir <- "../data/";
gf1   <- paste0(g.dir, "cancer_census_gene_tie1_2018.tsv");
gene1 <- read.table(gf1, header=T, sep="\t");
gf2   <- paste0(g.dir, "cancer_census_gene_tie2_2018.tsv");
gene2 <- read.table(gf2, header=T, sep="\t");
cgene <- unique(c(gene1[,1], gene2[,1]))

pid     <- "hsa04014";
pdat    <- read.table(system.file('extdata', 'relation_pathway_hsa.txt', package = "OmicPath"), header=T);
pdat.i  <- which(pdat[,3] == pid);
pdat    <- pdat[pdat.i,]; 
pdat.j  <- which(pdat[,1]==pdat[,2]);
if (length(pdat.j) > 0){
  pdat    <- pdat[-pdat.j,];
}

ct   <- "BRCA";
cols <- rainbow(10, alpha=0.9);
col1 <- rainbow(10, alpha=0.2);
col2 <- rainbow(10, alpha=0.8);

f.s  <- ct;
## cluster file
infile <- paste0("../data/clu_", f.s, ".txt.gz");
clu.d  <- read.table(gzfile(infile), header=T);
clu.u  <- unique(clu.d[,2]);
pdff   <- paste0("", f.s, ".pdf");
pdff3  <- paste0("", f.s, "_v3.pdf");
outf   <- paste0("", f.s, "_out.txt");
outf2  <- paste0("", f.s, "_net_out.txt");
out.s  <- NULL;
out.n  <- NULL;

## cluster
out.s <- NULL;
out.n <- NULL;
for (clu.i in clu.u){
  g.i <- which(clu.d[,2]==clu.i);
  gen <- clu.d[g.i,1];
  ## map genes on pathway
  out  <- pathway_map(gen, pdat);
  if (length(out)<2){
    next;
  }
  print(clu.i);
  tmp1 <- out[[1]];
  tmp2 <- out[[2]];
  df.s <- data.frame(rep(clu.i, nrow(tmp1)), tmp1);
  df.n <- data.frame(rep(clu.i, nrow(tmp2)), tmp2);
  out.s <- rbind(out.s, df.s);
  out.n <- rbind(out.n, df.n);
}
colnames(out.s) <- c("cluster", "sub", "gene");
colnames(out.n) <- c("cluster",	"sub", "from", "to", "relation.type");
write.table(out.s, outf, sep="\t", quote=F, row.names=F);
write.table(out.n, outf2, sep="\t", quote=F, row.names=F);

## plot network
pdf(pdff, 10, 5);
par(mfrow=c(1, 2));
clu.n <- paste0(out.n[,1], ".", out.n[,2]);
cl.u  <- unique(clu.n);
for (c.j in cl.u){
  out.j <- which(clu.n==c.j);
  out.x <- out.n[out.j,];
  print(nrow(out.x))
  if (nrow(out.x) < 3){
    next;
  }
  out.y <- out.x[,-c(1,2,6)];
  g     <- graph_from_data_frame(out.y, directed = TRUE, vertices = NULL);
  plot_network(g, cgene);
}
dev.off();

## plot pathway
pdf(pdff3, 12.8, 7.2);
par(mar=c(2, 2, 2, 2));
out1  <- KEGGplot(pid=pid, Taxonomy="hs", plot=T);
clu.s <- paste0(out.s[,1], ".", out.s[,2]);
c.u   <- unique(clu.s);
if (ct == "LUAD"){
  coln  <- rainbow(10, alpha=0.3)[c(1,2,7,9)]; 
} else {
  coln  <- rainbow(10, alpha=0.3)[c(1,9,7)]; 
}

k     <- 0;
for (c.i in c.u){
  k   <- k + 1;
  c.i <- which(clu.s==c.i);
  g.n <- out.s[c.i,3];
  c.j <- which(out1[,1] %in% g.n);
  points(as.numeric(out1[c.j,2]), as.numeric(out1[c.j,3]), pch=19, col=coln[k], cex=5);
}
dev.off();

