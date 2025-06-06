% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_coSet.R
\name{coSet}
\alias{coSet}
\title{coSet}
\usage{
coSet(sets, verbose = T, numCores = 1)
}
\arguments{
\item{sets}{a list containing series of gene symbols}
}
\value{
a count matrix
}
\description{
Estimate co-expression genes distribution
}
\examples{
set.seed(4090)
elements <- 1:6
sets <- list(
  sample(elements, 3),
  sample(elements, 4),
  sample(elements, 2),
  sample(elements, 5),
  sample(elements, 3))
resM <- coSet(sets)
print(resM)
}
\author{
Weibin Huang<\email{hwb2012@qq.com}>
}
