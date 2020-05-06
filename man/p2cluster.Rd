\name{p2cluster}
\alias{p2cluster}

\title{
 do cluster from data frame
}

\description{
 do cluster of data frame
}

\usage{
  p2cluster(dat=dat, method="louvain", wt=NULL){
}

\arguments{
   \item{dat}{data.frame or matrix: column 1 and 2 are the genes from association analysis. Column 3 is weight.}
   \item{method}{should be one of louvain, fast_greedy, infomap,
                  label_prop, spinglass and weight}
   \item{wt}{cutoff of weihgt}                
}

\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}

\references{
##
}

\examples{
##
}

