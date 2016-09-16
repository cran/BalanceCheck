\name{getD}
\alias{getD}
\title{Get intermediate results for CrossNN}
\description{This function calculates D12 and D21 for the CrossNN test.}
\usage{
getD(distM,treated.index)
}
\arguments{
  \item{distM}{The distance matrix for the pooled observations (pooled over the treated subjects and the matched controls). If there are n treated subjects and n matched controls, then this distance matrix is a 2n by 2n matrix with the [i,j] element the distance between observation i and observation j.  What distance to use is decided by users. Some simple choices are the Euclidean distance, L1 distance, and mahalanobis distance.}
  \item{treated.index}{The subject indices of the treated subjects.  The subjects are ordered in the same way as for calculating the distance matrix, distM.}
}
\seealso{
  \code{\link{CrossNN}}
}

