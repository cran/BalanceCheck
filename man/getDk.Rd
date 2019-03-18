\name{getDk}
\alias{getDk}
\title{Get intermediate results for CrossNN under k-NN}
\description{This function calculates D11 and D22 for the CrossNN test when k-NN (k>1) is used.}
\usage{
getDk(distM,treated.index,k)
}
\arguments{
  \item{distM}{The distance matrix for the pooled observations (pooled over the treated subjects and the matched controls). If there are n treated subjects and n matched controls, then this distance matrix is a 2n by 2n matrix with the [i,j] element the distance between observation i and observation j.  What distance to use is decided by users. Some simple choices are the Euclidean distance, L1 distance, and mahalanobis distance.}
  \item{treated.index}{The subject indices of the treated subjects.  The subjects are ordered in the same way as for calculating the distance matrix, distM.}
  \item{k}{An integer larger than 1.}
}
\seealso{
  \code{\link{CrossNN}}
}

