% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_GSClassifier_Data.R
\name{GSClassifier_Data}
\alias{GSClassifier_Data}
\title{Data in \code{GSClassifier} package}
\usage{
GSClassifier_Data(model = T)
}
\arguments{
\item{model}{logic. Whether to only show a list of \code{GSClassifier} models}
}
\description{
Internal dataset of GSClassifier package
}
\details{
ImmuneSubtype is a model for personalized calling of immune subtype by Thorsson et al. in 2018. \cr
PAD.train_20200110 is a model for personalized calling of Pan-immune Activation Dysfunction Subtypes(PAD).
}
\examples{
PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier"))
ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
