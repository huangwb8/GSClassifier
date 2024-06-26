% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelSeries_fitModels2.R
\name{fitEnsembleModel}
\alias{fitEnsembleModel}
\title{fitEnsembleModel}
\usage{
fitEnsembleModel(
  Xs,
  Ys,
  geneSet = NULL,
  na.fill.method = c("quantile", "rpart", NULL)[1],
  na.fill.seed = 2022,
  n = 20,
  sampSize = 0.7,
  sampSeed = 2020,
  breakVec = c(0, 0.25, 0.5, 0.75, 1),
  params = list(device = "cpu", nrounds = 15, max_depth = 10, eta = 0.5, nthread = 5,
    colsample_bytree = 1, min_child_weight = 1),
  nround.mode = c("fixed", "polling")[2],
  xgboost.seed = 105,
  caret.grid = expand.grid(nrounds = c(10, 15), max_depth = c(5, 10), eta = c(0.01, 0.1,
    0.3), gamma = c(0.5, 0.3), colsample_bytree = 1, min_child_weight = 1, subsample =
    0.7),
  caret.seed = 101,
  ptail = 0.5,
  verbose = F,
  numCores = 2
)
}
\arguments{
\item{Xs}{Gene expression matrix.}

\item{Ys}{Phenotype vector, multiclass}

\item{geneSet}{A list of genes for classification}

\item{na.fill.method}{Missing value imputation method for \code{\link{na_fill}} function. One of \code{'quantile'}, \code{'rpart'} and \code{NULL}.}

\item{na.fill.seed}{Seed for \code{\link{na_fill}} function.}

\item{n}{Size of the ensember, where each member is a result from
fitSubtypeModel}

\item{sampSize}{proportion of samples to hold back}

\item{sampSeed}{random seed for subset of Xs}

\item{breakVec}{vector of break points, used to bin expression data}

\item{params}{The parameters for \code{\link[xgboost]{xgb.train}}. 1. xgb.cv only: nfold; 2. xgboost: nrounds, max_depth, eta, nthread, colsample_bytree, min_child_weight.}

\item{nround.mode}{One of \code{fixed} and \code{polling}. \code{fixed} mode is recommended!
\itemize{
  \item \code{polling} Default but legacy feature, which means to call the \code{best_iteration} via \code{xgb.cv}
  \item \code{fixed} Use the default \code{nrounds} in \code{params}, so it's faster(10-20 folds) than \code{polling}.
}}

\item{xgboost.seed}{Seed for xgboost.}

\item{caret.seed}{The random seed for caret::train process when \code{params} is \code{NULL}}

\item{ptail}{Binary phenotype vector.}

\item{verbose}{whether report modeling process}
}
\value{
A list of lists of xgboost classifiers
}
\description{
Train multiple subtype models using cross validation
}
\details{
The geneid of \code{geneSet} and \code{Xs} must be the same (one of
  ENSEMBL, SYMBOL or ENTREZID). In addition, if \code{fitEnsembleModel} is hanged on, please check the space use of \code{/}, which would disturb the work of \code{\link[parallel]{makeCluster}}.
}
