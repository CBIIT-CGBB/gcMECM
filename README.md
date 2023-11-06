# gcMECM: Graph Clustering of Mutual Exclusivity of Cancer Mutations 

The package constructs the mutually exclusively mutated gene networks from mutation associations and gene interactions with the graph clustering technique 
and **identifies sub-networks with distinct biological functions with the canonical pathways**. The sub-networks reveal 
crucial genes in the canonical pathway related to cancers, which could used to build better prediction models for clinical response and survival. 
gcMECM also provides informative visualization functionality of mutual exclusivity and network. The associated publication was **"gcMECM: graph clustering of mutual exclusivity of cancer mutations"** on [BMC Bioinformatics 2021](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04505-w).


For additional visualization and clustering analysis, the packages, ([NCIRASPathway](https://github.com/CBIIT-CGR/NCIRASPathway), [OmicPath](https://github.com/CBIIT-CGR/OmicPath), [scCorr](https://github.com/CBIIT-CGR/scCorr), [SubPath](https://github.com/CBIIT-CGBB/SubPath), and [GCluster](https://github.com/CBIIT-CGR/GCluster)), could be used. 

# Introduction

### Mutation and pathway data
The mis-sense mutation and clinical outcome datasets for BRCA (Breast invasive carcinoma) in The Cancer Genome Atlas (TCGA) were obtained from [The NCI Genomic Data Commons (GDC, version 6)](https://portal.gdc.cancer.gov). RAS pathway v2.0 is obtained from [NCI Ras Initiative](https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2015/ras-pathway-v2). The pathway structure and gene coordinates were created manually for the visualization. KEGG pathway images and gene relationships were from [KEGG database](https://www.genome.jp/kegg/pathway.html).

### Detection of modules with mutually exclusive mutations 
#### Step 1
Generate the pairwise gene-gene adjacency distance matrix from the p-value of one-tailed Fisher’s exact test and generalized linear models (glm). Select gene pairs below the p-value cutoff and negative correlation in glm. The Fisher’s exact test p values are used as the distance in the matrix.
#### Step 2
Convert the distance matrix into a weighted graph or network using R package igraph. The resulting graph was clustered into modules with the Louvain algorithm. 
#### Step 3
Overlay the modules onto the canonical pathways.

# Tutorial and examples
### Installation
``` r
library(devtools)
install_github("CBIIT-CGBB/gcMECM")
``` 
  
### Cluster using Fisher's p values
```r
library(gcMECM)

## load the output of Fisher's test
dat    <- read.table(gzfile(infile), header=T);

## convert the p values to distance matrix
dat.d  <- p2dist(dat.s[,c(gene.i,gene.j)], as.numeric(dat.s[,pvalue.l]));

## construct clusters from the distance matrix 
clu    <- dist2cluster(dat.d, wt=wt, method="louvain");
```
Download the exmaple codes ([02p2cluster.R](examples/02p2cluster.R))
### Mapping clusters on the pathways
```r

library(gcMECM);
## load the pathway data
library(NCIRASPathway);

## retrieve genes in the pathway
g.xy  <- get_node_layout();

## load gene relationships in the pathway
pdat  <- get_relations()

## read cluster data from dist2cluster
clu.d  <- read.table(gzfile(infile), header=T);
```
Download the example codes ([R codes](examples/03cluster_map.R)) for NCI RAS pathway in [NCIRASPathway](https://github.com/CBIIT-CGR/NCIRASPathway) package. The output figures are as the follows (three sub-networks in the pathway and plots).

<img src="examples/03_1cluster_map.png" width="400" height="240">  <img src="examples/03_2cluster_map.png" width="90" height="240">
  
The option example codes ([R codes](examples/05kegg_pathway.R)) were for Ras signaling pathway of KEGG with [OmicPath](https://github.com/CBIIT-CGR/OmicPath) package. The output figures are as the follows (two sub-networks in the pathway and plots).

<img src="examples/05kegg_map.png" width="400" height="210">  <img src="examples/05kegg_network.png" width="120" height="210">

# Advanced (optional) steps
### Mutually exclusive mutation plots
For mutation plot ([R codes](examples/04_1plot_ME.R)) and sub-networks ([R codes](examples/04_2plot_ME.R))

<img src="examples/04_1plot_ME.png" width="350" height="250">  <img src="examples/04_2plot_ME.png" width="90" height="250">


