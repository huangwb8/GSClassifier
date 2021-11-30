

#' @title callEnsemble
#' @description Make subtype calls for each sample
#' @param X a numeric RNA expression vector or matrix with gene names
#' @inheritParams callSubtypes
#' @inheritParams geneMatch
#' @param subtype Default subtype methods. Now, only \code{'PAD'} and
#'   \code{'ImmuneSubtype'} are available. ATTENTION: If you use self-defined
#'   data (ens, geneAnnoation, geneSet or scaller), you MUST set \code{subtype}
#'   as \code{NULL}! All available options please visit
#'   \code{\link{GSClassifier_Data}}
#' @param scaller scaller object. A \code{\link[xgboost]{xgb.train}} object.
#' @import xgboost
#' @param ens A List of multiple model. Result of \code{\link{fitEnsembleModel}}
#' @return table, column 1 is best call, remaining columns are subtype
#'   prediction scores.
#' @examples
#' ## X is a matrix with SYMBOL ID in row scale.
#' resIS <- callEnsemble(
#' X,
#' ens = NULL,
#' geneAnnotation = NULL,
#' geneSet = NULL,
#' geneids = "symbol",
#' nClust = 6,
#' subtype = 'ImmuneSubtype',
#' bestMethod = c("maximum", "scaller")[2],
#' scaller = NULL
#' )
#' @export
callEnsemble <- function(X,
                         ens = NULL,
                         geneAnnotation = NULL,
                         geneSet = NULL,
                         scaller = NULL,
                         geneid = 'ensembl',
                         subtype = c('PAD.train_20200110', 'ImmuneSubtype')[1],
                         verbose = T) {
  X <- rightX(X)


  if (is.vector(X)) {
    # Single sample. Use one-type functions
    callEnsemble_One(
      X = X,
      ens = ens,
      geneAnnotation = geneAnnotation,
      geneSet = geneSet,
      scaller = scaller,
      geneid = geneid,
      subtype = subtype,
      verbose = verbose
    )
  } else {
    # Multiple samples. Use multi-type functions
    callEnsemble_Multi(
      X = X,
      ens = ens,
      geneAnnotation = geneAnnotation,
      geneSet = geneSet,
      scaller = scaller,
      geneid = geneid,
      subtype = subtype,
      verbose = verbose
    )

  }

}


#' @description Make subtype calls for single sample
#' @inheritParams callEnsemble
#' @import xgboost
#' @return table, column 1 is best call, remaining columns are subtype
#'   prediction scores.
callEnsemble_One <- function(X,
                             ens = NULL,
                             geneAnnotation = NULL,
                             geneSet = NULL,
                             scaller = NULL,
                             geneid = 'ensembl',
                             subtype = c('PAD.train_20200110', 'ImmuneSubtype')[1],
                             verbose = T) {
  ## Test
  if (F) {
    X = expr2[, 1]
    names(X) <- rownames(expr2)
    # X = expr2[,1:2]
    ens = NULL
    geneAnnotation = NULL
    geneSet = NULL
    scaller = NULL
    geneid = 'ensembl'
    subtype = 'PAD.train_20200110'
    verbose = T
  }

  ## Load model data
  if (!is.null(subtype)) {
    # Use system data
    if (verbose)
      LuckyVerbose('Use ', subtype, ' classifier...')
    l <-
      readRDS(system.file("extdata", paste0(subtype, '.rds'), package = "GSClassifier"))
    scaller <- l$scaller$Model
    geneAnnotation <- l$geneAnnotation
    ens = l$ens$Model
    nClust = length(ens[[1]])
    geneSet = l$geneSet
  } else {
    # Use self-defined data
    if (verbose)
      LuckyVerbose('Use self-defined classifier...')
    nClust = length(ens[[1]])
  }

  ## Matched data
  res0 <- geneMatch(X, geneAnnotation, geneid)
  X2 <- data.frame(target = res0$Subset, target2 = res0$Subset)
  matchError <- res0$matchError
  reportError(matchError)

  ## Call subtypes
  eList <-
    lapply(ens, function(ei)
      callSubtypes(mods = ei, X = X2, geneSet, nClust, verbose))
  ePart <- lapply(eList, function(a)
    a[, 3:(2 + nClust)])
  eStack <- array(unlist(ePart) , c(ncol(X2), nClust, length(ens)))
  eMeds  <- apply(eStack , 1:2 , median)
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:nClust # names(mods)

  ## Best call of maximum strategy
  bestCall_max <-
    apply(eMeds, 1, function(pi)
      colnames(eMeds)[which(pi == max(pi)[1])])

  ## Best call based on scaller
  if(!is.null(scaller)){
    bestCall_sc <- predict(scaller, as.matrix(eMeds)) + 1
  } else {
    bestCall_sc <- NA
  }


  ## Merge
  sampleIDs <- eList[[1]][, 1]
  res0 <- data.frame(
    SampleIDs = sampleIDs,
    BestCall = bestCall_sc,
    BestCall_Max = bestCall_max,
    eMeds
  )
  colnames(res0)[4:(3 + nClust)] <- 1:nClust
  LuckyVerbose('All done!')
  return(res0[1, ])
}


#' @description Make subtype calls for multiple samples
#' @inheritParams callEnsemble
#' @import xgboost
#' @return table, column 1 is best call, remaining columns are subtype
#'   prediction scores.
callEnsemble_Multi <- function(X,
                               ens = NULL,
                               geneAnnotation = NULL,
                               geneSet = NULL,
                               scaller = NULL,
                               geneid = 'ensembl',
                               subtype = c('PAD.train_20200110', 'ImmuneSubtype')[1],
                               verbose = T) {
  ## Load model data
  if (!is.null(subtype)) {
    # Use system data
    LuckyVerbose('Use ', subtype, ' classifier...')
    l <-
      readRDS(system.file("extdata", paste0(subtype, '.rds'), package = "GSClassifier"))
    scaller <- l$scaller$Model
    geneAnnotation <- l$geneAnnotation
    ens = l$ens$Model
    nClust = length(ens[[1]])
    geneSet = l$geneSet
  } else {
    # Use self-defined data
    LuckyVerbose('Use self-defined classifier...')
    nClust = length(ens[[1]])
  }

  ## Matched data
  res0 <- geneMatch(X, geneAnnotation, geneid)
  X <- res0$Subset
  matchError <- res0$matchError
  reportError(matchError)

  ## Call subtypes
  eList <-
    lapply(ens, function(ei)
      callSubtypes(mods = ei, X = X, geneSet, nClust, verbose))
  ePart <- lapply(eList, function(a)
    a[, 3:(2 + nClust)])
  eStack <- array(unlist(ePart) , c(ncol(X), nClust, length(ens)))
  eMeds  <- apply(eStack , 1:2 , median)
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:nClust # names(mods)

  ## Best call of maximum strategy
  bestCall_max <-
    apply(eMeds, 1, function(pi)
      colnames(eMeds)[which(pi == max(pi)[1])])

  ## Best call based on scaller
  if(!is.null(scaller)){
    bestCall_sc <- predict(scaller, as.matrix(eMeds)) + 1
  } else {
    bestCall_sc <- NA
  }

  ## Merge
  sampleIDs <- eList[[1]][, 1]
  res0 <- data.frame(
    SampleIDs = sampleIDs,
    BestCall = bestCall_sc,
    BestCall_Max = bestCall_max,
    eMeds
  )
  colnames(res0)[4:(3 + nClust)] <- 1:nClust
  LuckyVerbose('All done!')
  return(res0)
}


#' @title parCallEnsemble
#' @description Parallel version of callEnsemble
#' @inheritParams callEnsemble
#' @param numCores No. of CPU core
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @import foreach xgboost
#' @details Data of one sample was not supported. Please use \code{\link{callEnsemble}} instead.
#' @return table, column 1 is best call, remaining columns are subtype
#'   prediction scores.
#' @export
parCallEnsemble <- function(X,
                            ens = NULL,
                            geneAnnotation = NULL,
                            geneSet = NULL,
                            scaller = NULL,
                            geneids = 'ensembl',
                            subtype = c('PAD.train_20200110', 'ImmuneSubtype')[1],
                            verbose = T,
                            numCores = 2) {
  ## Load classifier data
  if (!is.null(subtype)) {
    # Use system data
    LuckyVerbose('Use ', subtype, ' classifier...')
    l <-
      readRDS(system.file("extdata", paste0(subtype, '.rds'), package = "GSClassifier"))
    scaller <- l$scaller$Model
    geneAnnotation <- l$geneAnnotation
    ens = l$ens$Model
    nClust = length(ens[[1]])
    geneSet = l$geneSet
  } else {
    # Use self-defined data
    LuckyVerbose('Use self-defined classifier...')
    nClust = length(ens[[1]])
  }

  ## Matched data and splited
  res0 <- geneMatch(X, geneAnnotation, geneids)
  X <- res0$Subset
  matchError <- res0$matchError
  reportError(matchError)
  XL <- spliteMatrix(X, cutoff = 10)

  ## Parallel call subtypes
  if (T) {
    # eList_bind
    eList_bind <- function(eL1, eL2) {
      if (is.null(eL1)) {
        return(eL2)
      } else {
        eL <- list()
        for (i in 1:length(eL1))
          eL[[i]] <- rbind(eL1[[i]], eL2[[i]])
        return(eL)
      }
    }

    # Parallel
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    clusterExport(cl, c('ens', 'XL', 'geneSet', 'nClust', 'verbose'), envir =
                    environment())
    clusterExport(
      cl,
      c(
        'callSubtypes',
        'callOneSubtype',
        'dataProc',
        'predict',
        'createPairsFeatures',
        'makeSetData',
        'breakBin',
        'binaryGene',
        'str_detect',
        'eList_bind'
      ),
      envir = environment()
    )
    eList <-
      foreach(i = 1:length(XL), .combine = eList_bind) %dopar% lapply(ens, function(ei)
        callSubtypes(mods = ei, X = XL[[i]], geneSet, nClust, verbose))
    stopImplicitCluster()
    stopCluster(cl)
  }

  ## Clean
  ePart <- lapply(eList, function(a)
    a[, 3:(2 + nClust)])
  eStack <- array(unlist(ePart) , c(ncol(X), nClust, length(ens)))
  eMeds  <- apply(eStack , 1:2 , median)
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:nClust # names(mods)

  ## Best call
  ## Best call of maximum strategy
  bestCall_max <-
    apply(eMeds, 1, function(pi)
      colnames(eMeds)[which(pi == max(pi)[1])])

  ## Best call based on scaller
  if(!is.null(scaller)){
    bestCall_sc <- predict(scaller, as.matrix(eMeds)) + 1
  } else {
    bestCall_sc <- NA
  }

  ## Merge
  sampleIDs <- eList[[1]][, 1]
  res0 <- data.frame(
    SampleIDs = sampleIDs,
    BestCall = bestCall_sc,
    BestCall_Max = bestCall_max,
    eMeds
  )
  colnames(res0)[4:(3 + nClust)] <- 1:nClust
  LuckyVerbose('All done!')
  return(res0)

}
