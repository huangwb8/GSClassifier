
####====================Basic function===================####

#' @title geneMatch
#' @description Match the incoming data to what was used in training
#' @param X gene expression matrix, genes in rows, samples in columns
#' @param geneAnnotation A data frame with ENSEMBL, SYMBOL and ENTREZID of the genes in \code{geneSet}
#' @param geneid type of gene id in \code{X}. One of \code{'symbol'}, \code{'entrez'} and \code{'ensembl'}.
#' @return Matched index of genes to the expression matrix.
#' @export
geneMatch <- function(X,
                      geneAnnotation = NULL,
                      geneid='ensembl') {

  ## Gene annotation
  if(is.null(geneAnnotation)){
    geneAnnotation <- readRDS(system.file("extdata", "modelGeneAnnotation.rds", package = "GSClassifier"))
  }

  name_col <- c('SYMBOL','ENSEMBL','ENTREZID'); names(name_col) <- c('symbol','ensembl','entrez')

  if (geneid == 'symbol') {
    idx <- match(table = rownames(X), x = as.character(geneAnnotation$SYMBOL))  ### this is just for the EBPP genes ###

  } else if (geneid == 'entrez') {
    idx <- match(table = rownames(X), x = as.character(geneAnnotation$ENTREZID))

  } else if (geneid == 'ensembl') {
    idx <- match(table = rownames(X), x = geneAnnotation$ENSEMBL)

  } else {
    print("For geneids, please use:  symbol, entrez, ensembl")
    return(NA)
  }

  matchError <- sum(is.na(idx)) / nrow(geneAnnotation)

  X2 <- X[idx,]  ### Adds NA rows in missing genes
  rownames(X2) <- as.character(geneAnnotation[,name_col[names(name_col) == geneid]])

  return(list(Subset=X2, matchError=matchError))
}


#' @description Error report
#' @param err the \code{matchError} result of \code{\link{geneMatch}}
reportError <- function(err) {
  print("**************************************")
  print("    Gene Match Error Report           ")
  print("                                      ")
  print(paste0("  percent missing genes: ",err*100,"           "))
  print("                                      ")
  print("**************************************")
}



#' @description given a matrix and gene pairs, create binary features
#' @param X One gene expression matrix
#' @param genes a vector of gene pairs
#' @return xbin, binned expression profile
#' @examples
#' Xsub <- createPairsFeatures(Xmat, genepairs)
createPairsFeatures <- function(X, genes) {

  # first convert the gene pairs to a named list
  # where each entry of the list is the genes for a given pivot-gene
  genePairs <- strsplit(genes, ':')
  pairList <- list()
  for(gi in genePairs) {
    pairList[[gi[1]]] <- c(pairList[[gi[1]]], gi[2])
  }

  resList <- list() # then for each pivot gene
  for (gi in names(pairList)) {
    # assuming it's in the data ... really should be!
    if (gi %in% rownames(X)) {
      gs <- pairList[[gi]]                  ## can end up with the pivot gene in the genes...
      pval <- as.numeric(X[gi,])            ## pivot values across samples
      idx <- match(table=rownames(X), x=gs) ## get index to genes for this pivot
      Xsub <- X[idx,]                       ## subset the matrix, NAs for missing genes, pivot gene on top
      if (class(Xsub) == 'numeric' & length(gs) == 1) {
        Xsub <- matrix(data=Xsub, ncol=ncol(X), nrow=1)
        colnames(Xsub) <- colnames(X)
      }
      rownames(Xsub) <- gs                  ## give gene IDs
      res0 <- lapply(1:ncol(Xsub), function(a) binaryGene(pval[a], Xsub[,a]))  ## create binary values
      resList[[gi]] <- do.call('cbind', res0)
    } else { # else we need to include some dummy rows
      randMat <- matrix(data=rbinom(n = length(pairList[[gi]]) * ncol(X), prob = 0.5, 1), ncol=ncol(X))
      colnames(randMat) <- colnames(X)
      resList[[gi]] <- randMat
    }
  }
  newMat <- do.call('rbind', resList)
  rownames(newMat) <- genes
  return(newMat)
}


#' @description splite matrix to a list of subset
#' @inheritParams geneMatch
#' @param cutoff the number of every subset of X
spliteMatrix <- function(X,cutoff=20){


  # Cut_vector
  cut_vector <- function(vt,nsplit=100){
    if(nsplit==1){
      v2 <- list(vt);
      names(v2) <- paste0("1-",length(v2))
    } else {
      ##转化成位置向量
      len.vt=1:length(vt)

      ##间隔
      len1 <- floor(length(len.vt)/nsplit)

      ## low.ci and upper.ci
      low.ci <- 1;
      for(i in 1:(nsplit-1)){
        low.ci[i+1] <- low.ci[i] + len1
      }
      upper.ci <- len1
      for(i in 1:(nsplit-1)){
        upper.ci[i+1] <- upper.ci[i] + len1
      }

      ## upper.ci的最后一个值
      upper.ci[length(upper.ci)] <- length(vt)

      ## 形成列表
      v2 <- NULL
      for(i in 1:length(low.ci)){
        v2.i <- low.ci[i]:upper.ci[i]
        v2.i <- vt[v2.i]
        v2 <- c(v2,list(v2.i))
        names(v2)[i] <- paste0(low.ci[i],"-",upper.ci[i])
      }
    }
    ## 输出结果
    return(v2)
  }

  # Splite matrix
  nXL <- floor(ncol(X)/cutoff)
  eXL <- cut_vector(colnames(X),nXL)
  XL <- list()
  for(i in 1:length(eXL)){ # i=1
    XL[[names(eXL)[i]]] <- X[,eXL[[i]]]
  }
  return(XL)
}



####==============callSubtype functions==================####

#' @description Data preprocessing
#' @param mods A model or list of models, containing breakpoints, used to bin expression data
#' @param geneSet A list of genes for classification
#' @param nClust The number of clustering
#' @inheritParams geneMatch
#' @importFrom stringr str_detect
#' @return Xbin, the binned, subset, and binarized values.
#' @examples
#' mod1 <- dataProc(X, mods)
dataProc <- function(X,
                     mods=NULL,
                     geneSet,
                     nClust=4) {

  # working with matrices
  Xmat <- as.matrix(X)

  # Set features
  nGS <- length(geneSet)
  featureNames <- c()
  for (j1 in 1:(nGS-1)) {
    for (j2 in (j1+1):nGS) {
      featureNames <- c(featureNames, paste0('s',j1,'s',j2))
    }
  }

  # get out the relevant items
  breakVec <- mods$breakVec
  genes    <- mods$bst$feature_names
  singleGenes <- genes[!str_detect(genes, ':')]
  singleGenes <- singleGenes[!singleGenes %in% featureNames]
  pairedGenes <- genes[str_detect(genes, ':')]

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec)
  rownames(Xbinned) <- rownames(Xmat)

  # and subset the genes to those not in pairs
  Xbinned <- Xbinned[singleGenes,]

  # here we have expression data, and we're using the pairs model
  # so we need to make pairs features.
  Xpairs <- createPairsFeatures(Xmat, pairedGenes)
  colnames(Xpairs) <- colnames(Xmat)

  # gene set features.
  Xset <- makeSetData(Xmat,geneSet)

  # join the data types and transpose
  Xbin <- t(rbind(Xbinned, Xpairs, Xset))
  return(Xbin)
}


#' @description Make subtype calls for one sample
#' @param mods xgboost model list
#' @param X gene expression matrix, genes in rows, samples in columns
#' @param ci cluster label, and index into mods
#' @inheritParams dataProc
#' @import xgboost
#' @return preds of one cluster model.
callOneSubtype <- function(mods, X, ci, geneSet, nClust) {

  # Xbin needs to have the same columns as the training matrix...
  print(paste0('calling subtype ', ci))
  mi <- mods[[ci]]
  Xbin <- dataProc(X, mods=mi, geneSet, nClust)
  pred <- predict(mi$bst, Xbin)
  return(pred)
}


#' @description Make subtype calls for each sample
#' @inheritParams callOneSubtype
callSubtypes <- function(mods, X, geneSet, nClust=4) {

  pList <- lapply(1:nClust, function(mi) callOneSubtype(mods, X, mi, geneSet, nClust))
  pMat  <- do.call('cbind', pList)
  colnames(pMat) <- 1:nClust # names(mods)
  bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])

  return(data.frame(SampleID=colnames(X), BestCall=bestCall, pMat, stringsAsFactors=F))
}


#' @title callEnsemble
#' @description Make subtype calls for each sample
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
callEnsemble <- function(
  X,
  ens=NULL,
  geneAnnotation = NULL,
  geneSet = NULL,
  scaller = NULL,
  geneid='ensembl',
  subtype = c('PAD.train_20200110',
              'PAD.all_20200110',
              'ImmuneSubtype')[1]
) {

  ## Load model data
  if(!is.null(subtype)){
    # Use system data
    cat('Use',subtype,'classifier...','\n')
    l <- readRDS(system.file("extdata", paste0(subtype,'.rds'), package = "GSClassifier"))
    scaller <- l$scaller$Model
    geneAnnotation <- l$geneAnnotation
    ens = l$ens$Model; nClust = length(ens[[1]])
    geneSet = l$geneSet
  } else {
    # Use self-defined data
    cat('Use self-defined classifier...','\n')
    nClust = length(ens[[1]])
  }

  ## Matched data
  res0 <- geneMatch(X,geneAnnotation,geneid)
  X <- res0$Subset
  matchError <- res0$matchError
  reportError(matchError)

  ## Call subtypes
  eList <- lapply(ens, function(ei) callSubtypes(mods=ei, X=X, geneSet, nClust))
  ePart <- lapply(eList, function(a) a[,3:(2+nClust)])
  eStack <- array( unlist(ePart) , c(ncol(X), nClust, length(ens)) )
  eMeds  <- apply( eStack , 1:2 , median )
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:nClust # names(mods)

  ## Best call
  bestCall_max <- apply(eMeds, 1, function(pi) colnames(eMeds)[which(pi == max(pi)[1])])
  bestCall_sc <- predict(scaller, as.matrix(eMeds)) + 1

  ## Merge
  sampleIDs <- eList[[1]][,1]
  res0 <- data.frame(SampleIDs=sampleIDs,
                     BestCall = bestCall_sc,
                     BestCall_Max=bestCall_max,
                     eMeds)
  colnames(res0)[4:(3+nClust)] <- 1:nClust
  return(res0)
}


#' @title parCallEnsemble
#' @description Parallel version: Make subtype calls for each sample
#' @inheritParams callEnsemble
#' @param numCores No. of CPU core
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @import foreach xgboost
#' @return table, column 1 is best call, remaining columns are subtype
#'   prediction scores.
#' @export
parCallEnsemble <- function(
  X,
  ens=NULL,
  geneAnnotation = NULL,
  geneSet = NULL,
  scaller = NULL,
  geneids='ensembl',
  subtype = c('PAD.train_20200110',
              'PAD.all_20200110',
              'ImmuneSubtype')[1],
  numCores=2) {

  ## Load classifier data
  if(!is.null(subtype)){
    # Use system data
    cat('Use',subtype,'classifier...','\n')
    l <- readRDS(system.file("extdata", paste0(subtype,'.rds'), package = "GSClassifier"))
    scaller <- l$scaller$Model
    geneAnnotation <- l$geneAnnotation
    ens = l$ens$Model; nClust = length(ens[[1]])
    geneSet = l$geneSet
  } else {
    # Use self-defined data
    cat('Use self-defined classifier...','\n')
    nClust = length(ens[[1]])
  }

  ## Matched data and splited
  res0 <- geneMatch(X,geneAnnotation,geneids)
  X <- res0$Subset
  matchError <- res0$matchError
  reportError(matchError)
  XL <- spliteMatrix(X,cutoff=10)

  ## Parallel call subtypes
  if(T){
    # eList_bind
    eList_bind <- function(eL1,eL2){
      if(is.null(eL1)){
        return(eL2)
      } else {
        eL <- list()
        for(i in 1:length(eL1)) eL[[i]] <- rbind(eL1[[i]],eL2[[i]])
        return(eL)
      }
    }

    # Parallel
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    clusterExport(cl, c('ens','XL','geneSet','nClust'),envir=environment())
    clusterExport(cl, c('callSubtypes', 'callOneSubtype', 'dataProc',  'predict', 'createPairsFeatures', 'makeSetData', 'breakBin', 'binaryGene', 'str_detect','eList_bind'),  envir=environment())
    eList <- foreach(i = 1:length(XL), .combine = eList_bind) %dopar% lapply(ens, function(ei) callSubtypes(mods=ei, X=XL[[i]], geneSet, nClust))
    stopImplicitCluster()
    stopCluster(cl)
  }

  ## Clean
  ePart <- lapply(eList, function(a) a[,3:(2+nClust)])
  eStack <- array( unlist(ePart) , c(ncol(X), nClust, length(ens)) )
  eMeds  <- apply( eStack , 1:2 , median )
  eMeds <- as.data.frame(eMeds)
  colnames(eMeds) <- 1:nClust # names(mods)

  ## Best call
  bestCall_max <- apply(eMeds, 1, function(pi) colnames(eMeds)[which(pi == max(pi)[1])])
  bestCall_sc <- predict(scaller, as.matrix(eMeds)) + 1

  ## Merge
  sampleIDs <- eList[[1]][,1]
  res0 <- data.frame(SampleIDs=sampleIDs,
                     BestCall = bestCall_sc,
                     BestCall_Max=bestCall_max,
                     eMeds)
  colnames(res0)[4:(3+nClust)] <- 1:nClust
  return(res0)

}

