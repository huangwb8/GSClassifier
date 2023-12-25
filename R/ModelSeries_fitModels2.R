


#' @description Train a single subtype model
#' @param Xs Gene expression matrix.
#' @param Ys Phenotype vector, multiclass
#' @param caret.seed The random seed for caret::train process when \code{params} is \code{NULL}
#' @param na.fill.method Missing value imputation method for \code{\link{na_fill}} function. One of \code{'quantile'}, \code{'rpart'} and \code{NULL}.
#' @param na.fill.seed Seed for \code{\link{na_fill}} function.
#' @param xgboost.seed Seed for xgboost.
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
fitSubtypeModel <- function(Xs,
                            Ys,
                            geneSet,
                            na.fill.method = c('quantile','rpart',NULL)[1],
                            na.fill.seed=2022,
                            breakVec = c(0, 0.25, 0.5, 0.75, 1.0),
                            params = list(
                              nrounds = 15,
                              max_depth = 10,
                              eta = 0.5,
                              nthread = 5,
                              colsample_bytree = 1,
                              min_child_weight = 1
                            ),
                            xgboost.seed = 105,
                            caret.grid = expand.grid(
                              nrounds = c(10,15),
                              max_depth = c(5,10),
                              eta = c(0.01, 0.1, 0.3),
                              gamma = c(0.5, 0.3),
                              colsample_bytree = 1,
                              min_child_weight = 1,
                              subsample = 0.7
                            ),
                            caret.seed = 101,
                            ptail = 0.05,
                            numCores = 2,
                            verbose = F) {
  modelList <- list()
  allLabels <- unique(Ys)

  # Set seeds
  set.seed(caret.seed)
  caret.seeds <- sample(1:100000, size = length(allLabels), replace = F)
  set.seed(xgboost.seed)
  xgboost.seeds <- sample(1:100000, size = length(allLabels), replace = F)

  # Missing value imputation
  Xs <- na_fill(Xs,
                method = na.fill.method,
                seed = na.fill.seed,
                verbose = verbose)

  # Model training
  for (i in 1:length(allLabels)) {
    # i=1

    yi = allLabels[i]
    if (verbose)
      LuckyVerbose(paste0('Subtype: ', yi, '  processing data...'))
    res0 <-
      trainDataProc(Xs,
                    Ys,
                    geneSet = geneSet,
                    subtype = yi,
                    ptail = ptail)
    dat  <- res0$dat
    Xbin <- dat$Xbin
    if (verbose)
      LuckyVerbose(paste0(
        'Training using ',
        ncol(Xbin),
        ' features and ',
        nrow(Xbin),
        ' samples...'
      ))

    if (!is.null(params)) {
      csfr <- cvFitOneModel(
        Xbin = dat$Xbin,
        Ybin = dat$Ybin,
        params = params,
        breakVec = breakVec,
        genes = dat$Genes,
        seed = xgboost.seeds[i],
        verbose = verbose
        )
    } else {
      # The caret::train strategy for params selection. Time consuming
      csfr <-
        cvFitOneModel2(
          dat$Xbin,
          dat$Ybin,
          breakVec,
          dat$Genes,
          seed = caret.seeds[i],
          caret.grid = caret.grid,
          verbose = verbose,
          numCores = numCores
        )
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
#' @importFrom parallel makeCluster clusterExport stopCluster parLapply
#' @details The geneid of \code{geneSet} and \code{Xs} must be the same (one of
#'   ENSEMBL, SYMBOL or ENTREZID). In addition, if \code{fitEnsembleModel} is hanged on, please check the space use of \code{/}, which would disturb the work of \code{\link[parallel]{makeCluster}}.
#' @return A list of lists of xgboost classifiers
#' @export
fitEnsembleModel <- function(Xs,
                             Ys,
                             geneSet = NULL,
                             na.fill.method = c('quantile','rpart',NULL)[1],
                             na.fill.seed=2022,
                             n = 20,
                             sampSize = 0.7,
                             sampSeed = 2020,
                             breakVec = c(0, 0.25, 0.5, 0.75, 1.0),
                             params = list(
                               nrounds = 15,
                               max_depth = 10,
                               eta = 0.5,
                               nthread = 5,
                               colsample_bytree = 1,
                               min_child_weight = 1
                             ),
                             xgboost.seed = 105,
                             caret.grid = expand.grid(
                               nrounds = c(10,15),
                               max_depth = c(5,10),
                               eta = c(0.01, 0.1, 0.3),
                               gamma = c(0.5, 0.3),
                               colsample_bytree = 1,
                               min_child_weight = 1,
                               subsample = 0.7
                             ),
                             caret.seed = 101,
                             ptail = 0.5,
                             verbose = F,
                             numCores = 2) {
  if (is.null(geneSet)) {
    geneSet = readRDS(system.file("extdata", paste0('PAD.train_20200110.rds'), package = "GSClassifier"))$geneSet
    LuckyVerbose('PAD subtype training...')
  }

  # Missing value imputation
  Xs <- na_fill(Xs,
                method = na.fill.method,
                seed = na.fill.seed,
                verbose = verbose)

  # Parallel Cores # https://stackoverflow.com/questions/21773199/r-makecluster-hangs-on-localhost-in-linux
  if (verbose) LuckyVerbose('Start parallel process...')
  cl <- makeCluster(numCores)

  # Seeds
  set.seed(sampSeed)
  seeds <- sample(1:100000, size = n, replace = F)
  set.seed(xgboost.seed)
  xgboost.seeds <- sample(1:100000, size = n, replace = F)

  # Assistant function
  fitFun <- function(i, verbose) {
    # i = (1:n)[1]
    modi <- c()
    set.seed(seeds[i])
    jdx <- sample(1:ncol(Xs), size = sampSize * ncol(Xs), replace = F)
    Xs2 <- Xs[, jdx]
    Ys2 <- Ys[jdx]
    modi <- fitSubtypeModel(
      Xs = Xs2,
      Ys = Ys2,
      geneSet = geneSet,
      na.fill.method = NULL,
      na.fill.seed = na.fill.seed,
      breakVec = breakVec,
      params = params,
      xgboost.seed = xgboost.seeds[i],
      caret.grid = caret.grid,
      caret.seed = caret.seed,
      ptail = ptail,
      verbose = verbose,
      numCores = numCores
    )
    return(modi)
  }

  # Parallel parameters
  clusterExport(cl, 'Xs',  envir = environment())
  clusterExport(cl, 'Ys',  envir = environment())
  clusterExport(cl, 'geneSet',  envir = environment())
  clusterExport(cl, 'na.fill.method',  envir = environment())
  clusterExport(cl, 'na.fill.seed',  envir = environment())
  clusterExport(cl, 'sampSize',  envir = environment())
  clusterExport(cl, 'seeds',  envir = environment())
  clusterExport(cl, 'caret.seed',  envir = environment())
  clusterExport(cl, 'breakVec',  envir = environment())
  clusterExport(cl, 'params',  envir = environment())
  clusterExport(cl, 'xgboost.seeds',  envir = environment())
  clusterExport(cl, 'ptail',  envir = environment())
  clusterExport(cl, 'verbose',  envir = environment())
  clusterExport(cl, 'caret.grid',  envir = environment())
  clusterExport(cl, 'numCores',  envir = environment())
  clusterExport(
    cl,
    c(
      'fitSubtypeModel',
      'trainDataProc',
      'cvFitOneModel',
      'cvFitOneModel2',
      'makeSetData',
      'makeGenePairs',
      'breakBin',
      'binaryGene',
      'featureSelection',
      'testFun'
    ),
    envir = environment()
  )
  clusterExport(cl, c('na_fill', 'is.one.na', 'rpart', 'predict'), envir =
                  environment())
  clusterExport(
    cl,
    c(
      'xgb.DMatrix',
      'xgb.cv',
      'xgboost',
      'trainControl',
      'train',
      'xgb.train',
      'LuckyVerbose'
    ),
    envir = environment()
  )
  clusterExport(cl, 'fitFun',  envir = environment())

  # Parallel process
  ens <- parLapply(
    cl = cl,
    X = 1:n,
    fun = function(x)
      fitFun(x, verbose = verbose)
  )

  # Output results
  res <- list(
    Repeat = list(
      Xs = Xs,
      Ys = Ys,
      geneSet = geneSet,
      n = n,
      sampSize = sampSize,
      sampSeed = sampSeed,
      breakVec = breakVec,
      params = params,
      caret.grid = caret.grid,
      caret.seed = caret.seed,
      ptail = ptail,
      numCores = numCores
    ),
    Model = ens
  )

  if (verbose)
    LuckyVerbose('End parallel process!')
  stopCluster(cl)

  return(res)
}
