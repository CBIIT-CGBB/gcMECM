\name{combine_test}
\alias{combine_test}

\title{
 combine p values
}

\description{
 combine p values
}

\usage{
  combine_test(p, weight, method = c("fisher", "z.transform", "logit"), hetero = FALSE, na.rm = FALSE)
}

\arguments{
   \item{p}{a vector of p values}
   \item{weight}{weights for p values as a vector. The length should be the same as the p values.}
   \item{method}{one of "fisher", "z.transform" and "logit"}
   \item{hetero}{TRUE is the heterogeneity should be taken into account, FALSE otherwise}
   \item{na.rm}{TRUE if the missing values should be removed from the data, FALSE otherwise}
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

