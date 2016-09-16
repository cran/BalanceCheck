\name{getR1R2}
\alias{getR1R2}
\title{Get intermediate results for CrossMST}
\description{This function calculates R1 and R2 for the CrossMST test.}
\usage{
getR1R2(E,treated.index)
}
\arguments{
  \item{E}{A matrix with the number of rows the number of edges in the MST and 2 columns.  Each row records the subject indices of the two ends of an edge in the MST.}
  \item{treated.index}{The subject indices of the treated subjects.}
}
\seealso{
  \code{\link{CrossMST}}
}

