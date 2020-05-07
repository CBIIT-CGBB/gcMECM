# gcMECM: Graphical Clustering of Mutual Exclusivity of Cancer Mutations to identify sub-biological functions from pathways

The package combines graph clustering, mutation association, and gene interaction 
by the mutually exclusively mutated gene sub-networks to identify sub-biological functions from 
pathways. The sub-networks reveal 
crucial genes in the canonical pathway and discover new cancer-relevant genes and relative biological 
functions, which could used to build better prediction models of clinical response and survival. 
gcMECM also provides informative visualization functionality of mutual exclusivity and network.

# Introduction

### Mutation and pathway data
The mis-sense mutation and clinical outcome datasets for LUAD (Lung Adenocarcinoma) and BRCA (Breast invasive carcinoma) in The Cancer Genome Atlas (TCGA) were obtained from The NCI Genomic Data Commons (GDC) (https://portal.gdc.cancer.gov, version 6). RAS pathway v2.0 is obtained from NCI Ras Initiative  (https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2015/ras-pathway-v2). The pathway structure and gene coordinates were created manually for the visualization. KEGG pathway images and gene relationships were from KEGG database (https://www.genome.jp/kegg/pathway.html).

### Detection of modules with mutually exclusive mutations 

Graph-based clustering algorthim was used to construct modules with mutually exclusive mutation patterns in the R platform. The process consists of three steps: First, generate the pairwise gene-gene adjacency distance matrix from the p-value of one-tailed Fisher’s exact test and generalized linear models (glm). Select gene pairs with p-value < 0.0001 and negative correlation in glm. The Fisher’s exact test p values are used as the distance in the matrix. Next, The adjacency distance matrix was converted into a weighted graph or network using R package igraph. The resulting graph was clustered into modules with the Louvain algorithm, which is an efficient graph-clustering method based on the modularity measure and a heuristic approach (Blondel, 2008), by applying the cutoff of edge values less than 1.00E-08 and 1.00E-12 in BRCA and LUAD respectively. In each module, genes are closely related to each other with negative relationships due to the mutually exclusive mutations. Modules with less that three genes are removed from further analysis. Lastly, each module was overlayed onto canonical pathways from NCI Ras pathway (https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2015/ras-pathway-v2) or KEGG to identify sub-networks, which composes of genes directly connected in both mutually exclusive modules and canonical pathways. A combined p-value for sub-networks were calculated with gene-pair fisher’s tests using combine.test function in R package survcomp. All sub-networks in the context of pathways were visualized with R package igraph and further examined for enrichment analysis with g:Profiler (https://biit.cs.ut.ee/gprofiler/gost) against Gene Ontology and KEGG databases and survival analysis with R survival package.

# Tutorial and examples
### Installation
``` r
library(devtools)
remotes::install_github("CBIIT-CGR/gcMECM")
``` 
  
### Cluster using Fisher's p values
```r
library(gcMECM)

## load the output of Fisher's test
dat    <- read.table(gzfile(infile), header=T);

## convert the p values to distance matrix
dat.d  <- p2dist(dat.s[,c(gene.i,gene.j)], as.numeric(dat.s[,pvalue.l]));

## cluster from the distance matrix 
clu    <- dist2cluster(dat.d, wt=wt, method="louvain");
```
[example code](R/p2dist.R)
### Mapping clusters on the pathways
```r

```

# Advanced (optional) steps

# References
