

#' @description Train a single subtype model
#' @param Xs Gene expression matrix.
#' @param Ys Phenotype vector, multiclass
#' @param caret.seed The random seed for caret::train process when \code{params} is \code{NULL}
#' @param na.fill NA filling. One of 'rpart' or \code{NULL}.
#' @inheritParams cvFitOneModel
#' @inheritParams trainDataProc
#' @return A list of xgboost classifiers, one for each subtype.
#' @examples
#' params=list(
#' max_depth = 2,
#' eta = 0.5,
#' nrounds = 100,
#' nthread = 5,
#' nfold=5)
#' @export
fitSubtypeModel <- function(Xs, Ys,
                            geneSet,
                            na.fill = c('rpart', NULL)[1],
                            breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                            params=list(max_depth = 2,
                                        eta = 0.5,
                                        nrounds = 100,
                                        nthread = 5,
                                        nfold=5),
                            caret.grid = NULL, caret.seed = 101,
                            ptail = 0.05,
                            verbose = F) {

  modelList <- list()
  allLabels <- unique(Ys)

  # Set seeds
  set.seed(caret.seed); caret.seeds <- sample(1:100000,size= length(unique(Ys)),replace = F)

  # NA filling with recursive partitioning and regression trees
  if(!is.null(na.fill)){
    # Do NA filling
    if(ncol(Xs)>=3 & na.fill == 'rpart'){
      if(verbose) LuckyVerbose('Missing data filled with RPART algorithm...')
      Xs <- na_fill(Xs, method="anova", na.action = na.rpart)
    }
  }


  # Model training
  for (i in 1:length(allLabels)) { # i=1

    yi = allLabels[i]
    if(verbose) LuckyVerbose(paste0('Subtype: ',yi, '  processing data...'))
    res0 <- trainDataProc(Xs, Ys, geneSet = geneSet,subtype=yi, ptail=ptail)
    dat  <- res0$dat; Xbin <- dat$Xbin
    if(verbose) LuckyVerbose(paste0('Training using ',ncol(Xbin),' features and ',nrow(Xbin),' samples...'))

    if(!is.null(params)){
      # Best iteration based on one available model
      csfr <- cvFitOneModel(dat$Xbin, dat$Ybin, params, breakVec, dat$Genes, verbose=verbose)
    } else {
      # The caret::train strategy for params selection. Time consuming
      csfr <- cvFitOneModel2(dat$Xbin, dat$Ybin, breakVec, dat$Genes, seed = caret.seeds[i], caret.grid = caret.grid, verbose=verbose)
    }

    modelList[[yi]] <- csfr
  }

  names(modelList) <- allLabels
  return(modelList)
}

#' @title fitEnsembleModel
#' @description Train multiple subtype models using cross validation
#' @inheritParams fitSubtypeModel
#' @param n Size of the ensember, where each member is a result from
#'   fitSubtypeModel
#' @param sampSize proportion of samples to hold back
#' @param sampSeed random seed for subset of Xs
#' @param numCores number of cores to use, one per ensemble member
#' @importFrom parallel makeCluster clusterExport stopCluster parLapply
#' @details The geneid of \code{geneSet} and \code{Xs} must be the same (one of
#'   ENSEMBL, SYMBOL or ENTREZID).
#' @return A list of lists of xgboost classifiers
#' @export
fitEnsembleModel <- function(Xs, Ys,
                             geneSet = NULL,
                             na.fill = c('rpart', NULL)[1],
                             n = 20,
                             sampSize = 0.7, sampSeed = 2020,
                             breakVec =c(0, 0.25, 0.5, 0.75, 1.0),
                             params = list(max_depth = 5,
                                           eta = 0.5,
                                           nrounds = 100,
                                           nthread = 5,
                                           nfold = 5),
                             caret.grid = NULL, caret.seed = 101,
                             ptail = 0.5,
                             verbose = F,
                             numCores = 2) {

  if(is.null(geneSet)){
    geneSet = readRDS(system.file("extdata", paste0('PAD.train_20200110.rds'), package = "GSClassifier"))$geneSet
    LuckyVerbose('PAD subtype training...')
  }

  # NA filling with recursive partitioning and regression trees
  if(!is.null(na.fill)){
    # Do NA filling
    if(ncol(Xs)>=3 & na.fill == 'rpart'){
      if(verbose) LuckyVerbose('Missing data filled with RPART algorithm...')
      Xs <- na_fill(Xs, method="anova", na.action = na.rpart)
    }
  }

  # Parallel Cores
  cl <- makeCluster(numCores,  outfile='')

  # Seeds
  set.seed(sampSeed); seeds <- sample(1:100000,size=n,replace = F)

  # Assistant function
  fitFun <- function(i, verbose) { # i = (1:n)[1]
    modi <- c()
    set.seed(seeds[i]); jdx <- sample(1:ncol(Xs), size = sampSize * ncol(Xs), replace=F)
    Xs2 <- Xs[,jdx]
    Ys2 <- Ys[jdx]
    modi <- fitSubtypeModel(
      Xs=Xs2, Ys=Ys2,
      geneSet = geneSet, na.fill = NULL,
      breakVec = breakVec,
      params = params,
      caret.grid = NULL, caret.seed = caret.seed,
      ptail= ptail,
      verbose = verbose
    )
    return(modi)
  }

  # Parallel parameters
  clusterExport(cl, 'Xs',  envir=environment())
  clusterExport(cl, 'Ys',  envir=environment())
  clusterExport(cl, 'geneSet',  envir=environment())
  clusterExport(cl, 'na.fill',  envir=environment())
  clusterExport(cl, 'sampSize',  envir=environment())
  clusterExport(cl, 'seeds',  envir=environment())
  clusterExport(cl, 'caret.seed',  envir=environment())
  clusterExport(cl, 'breakVec',  envir=environment())
  clusterExport(cl, 'params',  envir=environment())
  clusterExport(cl, 'ptail',  envir=environment())
  clusterExport(cl, 'verbose',  envir=environment())
  clusterExport(cl, 'caret.grid',  envir=environment())
  clusterExport(cl, c('fitSubtypeModel','trainDataProc','cvFitOneModel','cvFitOneModel2','makeSetData','makeGenePairs','breakBin','binaryGene','featureSelection','testFun'),envir=environment())
  clusterExport(cl, c('na_fill','is.one.na','rpart','predict'),envir=environment())
  clusterExport(cl, c('xgb.DMatrix','xgb.cv','xgboost','trainControl','train','xgb.train','LuckyVerbose'),envir=environment())
  clusterExport(cl, 'fitFun',  envir=environment())

  # Parallel process
  if(verbose) LuckyVerbose('Start parallel process...')
  ens <- parLapply(cl=cl, X=1:n,
                   fun = function(x)fitFun(x,verbose=verbose))

  # Output results
  res <- list(
    Repeat = list(
      Xs = Xs, Ys = Ys,
      geneSet = geneSet,
      n=n,
      sampSize=sampSize, sampSeed = sampSeed,
      breakVec=breakVec,
      params=params,
      caret.grid = caret.grid,
      caret.seed = caret.seed,
      ptail=ptail,
      numCores=numCores
    ),
    Model = ens
  )

  if(verbose) LuckyVerbose('End parallel process!')
  stopCluster(cl)

  return(res)
}
