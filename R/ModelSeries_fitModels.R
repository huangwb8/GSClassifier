


####====================Basic function==================####

#' @description Get difference in mean rank sums (Ybin=0 vs. 1) for a single
#'   gene
#' @param rankg Gene expression profile for single sample
#' @param Ybin Binary phenotype vector.
#' @return test result, numeric value, the rank sum test.
#' @details The larger the absolute difference across phenotype, the more
#'   defferential distribution of the gene, which indicates that the gene could
#'   be a potential marker for prediciton of the phenotype. This is why
#'   \code{ptail} is introduced in \code{\link{featureSelection}} for
#'   improvement of calculation.
#' @examples
#' res1 <- testFun(G, Ybin)
#' @export
testFun <- function(rankg, Ybin) {
  res0 <- (sum(rankg[Ybin == 0], na.rm = T) / sum(Ybin == 0, na.rm = T)) - (sum(rankg[Ybin == 1], na.rm = T) / sum(Ybin == 1, na.rm = T))
  return(res0)
}


#' @description Subset the genes, given a matrix
#' @export
#' @param Xmat Matrix of gene expression data.
#' @param Ybin Binary phenotype vector.
#' @param testRes Result of \code{\link{testFun}} function.
#' @param ptail Proportion tail. Range: 0< ptail <=0.5。A larger \code{ptail}
#'   means longer process time but higher accurate prediction.
#' @return Xsub, subset of Xmat by genes
#' @examples
#' Xsub <- featureSelection(Xmat, Ybin, 0.1)
featureSelection <- function(Xmat, Ybin,
                             testRes,
                             ptail=0.5) {
  idx <- which( (testRes < quantile(testRes, ptail, na.rm = T)) |
                  (testRes > quantile(testRes, 1.0-ptail, na.rm = T)) )
  Xsub <- Xmat[idx,]
  Xsub[is.na(Xsub)] <- 0 # NA value would be turned as 0
  Xgenes <- rownames(Xmat)[idx]
  return(list(Xsub=Xsub, Genes=Xgenes))
}


#' @description Subset the genes, given a matrix
#' @param x One gene expression sample profile.
#' @param breakVec A vector of probabilities
#' @return xbin, binned expression profile
#' @examples
#' Xsub <- featureSelection(Xmat, Ybin, 0.1)
#'
breakBin <- function(x, breakVec){
  brks <- quantile(as.numeric(x), probs=breakVec, na.rm = T)
  xbin <- .bincode(x = x, breaks = brks, include.lowest = T)
  xbin <- as.numeric(xbin)
  xbin
}


#' @description Assistant function of \code{\link{makeGenePairs}} function
#' @details Give a vector and pivotvalue. If the element in the vetor is smaller
#'   than the pivotvalue, output \code{1}, else output \code{0}. Then, replace
#'   NAs with random values. This is why NA value in the expression matrix is
#'   agreeable.
#' @examples
#' x <- c(4.318257,9.456551,10.607902,6.716920,5.249954,8.482040,6.551659,8.912432,5.261267,3.624650,7.155214,8.608159,4.012154,6.344634,5.454208,7.119977,4.144586,5.727500,9.828227,9.692067,6.729285,8.426716,9.852803,10.267834,7.254560,11.441216,11.206160,10.838834,9.934967,10.079822,6.154237)
#' binaryGene(7.06627,x) # 1 0 0 1 1 0 1 0 1 1 0 0 1 1 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 1
binaryGene <- function(pivotvalue, values) {
  res0 <- sapply(values, function(b) as.numeric(b <= pivotvalue))
  res0[is.na(res0)] <- rbinom(n = sum(is.na(res0)), prob = 0.5, 1)  ## replace NAs with random values
  return(res0)
}


#' @description Make gene pairs. Assistant function.
#' @param genes Genes for pair
#' @param Xsub Subset of X
#' @details For every g1:g2 pair across samples, res_g1-g2 = expression of g1 -
#'   expresion of g2 was calculated. If res_g1-21 >0, output \code{1}, else
#'   \code{0}.
makeGenePairs <- function(genes, Xsub) {
  # collect results here
  resList <- list()
  # all pairs of genes here
  all.combos <- t(combn(genes,2))

  for (i in 1:nrow(all.combos)) {
    gi <- all.combos[i,1]  # the pivot gene
    gval <- as.numeric(Xsub[gi,]) # and it's value
    # get the paired genes out.
    gcom <- all.combos[all.combos[,1] == gi,2]
    # pick sample j, get pivot value j which is for gene i, and pivot all the paired genes
    res0 <- lapply(1:ncol(Xsub), function(j) binaryGene(gval[j], Xsub[gcom,j]))
    # make matrix of features.
    resMat <- do.call('cbind', res0)
    rownames(resMat) <- sapply(gcom, function(gj) paste0(gi,':',gj))
    resList[[gi]] <- resMat
  }
  Xbin <- do.call('rbind', resList)
  colnames(Xbin) <- colnames(Xsub)
  return(Xbin)
}


#' @description Make set data.
#' @param Xmat Matrix of gene expression data.
#' @param geneSet List of genes for classification.
makeSetData <- function(Xmat,geneSet) {

  resultList <- list()
  nGS <- length(geneSet)

  featureNames <- c()
  for (j1 in 1:(nGS-1)) {
    for (j2 in (j1+1):nGS) {
      featureNames <- c(featureNames, paste0('s',j1,'s',j2))
    }
  }

  # for each sample
  for (i in 1:ncol(Xmat)) {
    res0 <- numeric(length=length(featureNames))
    idx <- 1
    for (j1 in 1:(nGS-1)) {
      for (j2 in (j1+1):nGS) {
        set1 <- geneSet[[j1]]
        set2 <- geneSet[[j2]]
        vals1 <- Xmat[rownames(Xmat) %in% set1,i]
        vals2 <- Xmat[rownames(Xmat) %in% set2,i]
        res1 <- sapply(vals1, function(v1) sum(v1 > vals2, na.rm=T))
        res0[idx] <- sum(res1, na.rm = T) / (length(vals1) * length(vals2))
        idx <- idx+1
      }
    }
    resultList[[i]] <- as.numeric(res0)
  }
  resMat <- do.call(cbind, resultList)
  colnames(resMat) <- colnames(Xmat)
  rownames(resMat) <- featureNames
  return(resMat)
}

####========================Fit models==================####

#' @title trainDataProc
#' @description Binary data preprocessing based on expression matrix and real clustering
#' @param Xmat Matrix of gene expression, samples in columns, genes in rows
#' @param Yvec Vector of phenotype, strings, 1, 2, etc
#' @inheritParams dataProc
#' @param subtype cluster name, string '1', as in Yvec
#' @param ptail Binary phenotype vector.
#' @param breakVec vector of break points, used to bin expression data
#' @return List of Xbin and Ybin, the binned, subset, and binarized values.
#' @export
trainDataProc <- function(Xmat, Yvec,
                          geneSet, subtype=1,
                          ptail=0.5,
                          breakVec=c(0, 0.25, 0.5, 0.75, 1.0)) {

  # Create the binary subtype identifier
  Ybin <- ifelse(Yvec == subtype, yes = 1, no=0)

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec) # bin each column
  rownames(Xbinned) <- rownames(Xmat)

  # rank the data, and use it for feature selection
  Xrank <- apply(Xmat, 2, rank)
  testRes <- sapply(1:nrow(Xrank), function(gi) testFun(as.numeric(Xrank[gi,]), Ybin))  # get genes with big rank diffs.

  # subset the expression data for pairs
  Xfeat <- featureSelection(Xmat, Ybin, testRes, ptail)  # subset genes
  Xpairs <- makeGenePairs(Xfeat$Genes, Xfeat$Xsub)

  # subset the binned genes
  Xbinned <- Xbinned[Xfeat$Genes,]

  # gene set features.
  Xset <- makeSetData(Xmat,geneSet)

  # join the data types and transpose
  Xbin <- t(rbind(Xbinned, Xpairs, Xset))
  genes <- colnames(Xbin)

  return(list(dat=list(Xbin=Xbin,
                       Ybin=Ybin,
                       Genes=genes),
              testRes=testRes,
              breakVec=breakVec))
}


#' @description Train a single subtype model using cross validation
#' @param Xbin Binned and filtered gene expression matrix.
#' @param Ybin Binned phenotype vector.
#' @param params The parameters for \code{\link[xgboost]{xgb.train}}
#' @param genes Genes for modeling
#' @inheritParams trainDataProc
#' @importFrom xgboost xgboost xgb.cv xgb.DMatrix
#' @return A single xgboost classifier.
cvFitOneModel <- function(Xbin, Ybin,
                          params=list(max_depth = 2,
                                      eta = 0.5,
                                      nrounds = 100,
                                      nthread = 5, nfold=5),
                          breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                          genes){
  dtrain <- xgb.DMatrix(Xbin, label = Ybin)

  # xgb.cv
  for(i in 1:10000){
    x <- tryCatch(
      cvRes <- xgb.cv(data = dtrain,
                      nrounds=params$nrounds,
                      nthread=params$nthread,
                      nfold=params$nfold,
                      max_depth=params$max_depth,
                      eta=params$eta,
                      early_stopping_rounds=2,
                      metrics = list("error", "auc"),
                      objective = "binary:logistic"),
      error = function(e)e)
    if('message' %in% names(x)){
      cat('Attention! AUC: the dataset only contains pos or neg samples. Repeat xgb.cv','\n')
      x_error <- x
    } else {
      cvRes <- x
      break
    }
  }
  cat('Best interation: ',cvRes$best_iteration,'\n')

  # xgboost via best interation
  bst <- xgboost(data = Xbin,
                 label = Ybin,
                 max_depth=params$max_depth,
                 eta=params$eta,
                 nrounds = cvRes$best_iteration,
                 nthread=params$nthread,
                 objective = "binary:logistic")

  return(list(bst=bst, breakVec=breakVec, genes=genes))
}

#' @importFrom xgboost xgb.train xgb.DMatrix
#' @importFrom caret trainControl train
cvFitOneModel2 <- function(Xbin, Ybin,
                           breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                           genes,
                           seed = 101){

  # Data
  x = Xbin
  y = factor(Ybin,levels = c(0,1))

  # Parameters of caret::train
  grid <- expand.grid(
    nrounds = c(10,15,20,30,40),
    colsample_bytree = 1,
    min_child_weight = 1,
    eta = c(0.01, 0.1, 0.3),
    gamma = c(0.5, 0.3),
    subsample = 0.7,
    max_depth = c(2,3)
  )

  cntrl <- trainControl(
    method = "cv",
    number = 5,
    verboseIter = TRUE,
    returnData = FALSE,
    returnResamp = "final"
  )

  ## caret::train process
  set.seed(seed)
  train.xgb <- train(
    x = x,
    y = y,
    trControl = cntrl,
    tuneGrid = grid,
    method = "xgbTree"
  )

  ## Optimized parameters
  train.mat <- xgb.DMatrix(data = x, label = Ybin)
  param <- as.list(train.xgb$bestTune)
  nrounds <- param$nrounds
  param2 <- c(list(objective = "binary:logistic",
                   booster = "gbtree",
                   eval_metric = "error"),
              param[-1])
  xgb.fit <- xgb.train(params = param2,
                       data = train.mat,
                       nrounds = param$nrounds)

  return(list(bst=xgb.fit, breakVec=breakVec, genes=genes))


}

#' @description Train a single subtype model using cross validation
#' @param Xs Gene expression matrix.
#' @param Ys Phenotype vector, multiclass
#' @param caret.seed The random seed for caret::train process when \code{params} is \code{NULL}
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
                            breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                            params=list(max_depth = 2,
                                        eta = 0.5,
                                        nrounds = 100,
                                        nthread = 5,
                                        nfold=5),
                            caret.seed = 101,
                            ptail=0.05) {

  modelList <- list()
  allLabels <- unique(Ys)

  set.seed(caret.seed); caret.seeds <- sample(1:100000,size= length(unique(Ys)),replace = F)


  for (i in 1:length( allLabels)) { # i=1

    yi = allLabels[i]
    print(paste0('Subtype: ',yi, '  processing data...'))
    res0 <- trainDataProc(Xs, Ys, geneSet = geneSet,subtype=yi, ptail=ptail)
    print(paste0('   training using ', dim(res0$dat$Xbin), ' features x samples'))
    dat  <- res0$dat
    if(!is.null(params)){
      csfr <- cvFitOneModel(dat$Xbin, dat$Ybin, params, breakVec, dat$Genes)
    } else {
      csfr <- cvFitOneModel2(dat$Xbin, dat$Ybin, breakVec, dat$Genes,seed = caret.seeds[i])
    }

    modelList[[yi]] <- csfr
  }

  names(modelList) <- allLabels
  return(modelList)
}

#' @title fitEnsembleModel
#' @description Train a single subtype model using cross validation
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
                             n=20,
                             sampSize=0.7, sampSeed = 2020,
                             breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                             params=list(max_depth = 5,
                                         eta = 0.5,
                                         nrounds = 100,
                                         nthread = 5,
                                         nfold=5),
                             caret.seed = 101,
                             ptail=0.5,
                             numCores=2) {

  if(is.null(geneSet)){
    geneSet = readRDS(system.file("extdata", paste0('PAD.train_20200110.rds'), package = "GSClassifier"))$geneSet
    cat('PAD subtype training...', '\n')
  }

  cl <- makeCluster(numCores,  outfile='')

  set.seed(sampSeed); seeds <- sample(1:100000,size=n,replace = F)

  fitFun <- function(i) { # i = (1:20)[1]
    modi <- c()
    # set.seed(seeds[1]); jdx <- sample(1:ncol(Xs), size = sampSize * ncol(Xs), replace=F)
    set.seed(seeds[i]); jdx <- sample(1:ncol(Xs), size = sampSize * ncol(Xs), replace=F)
    Xs2 <- Xs[,jdx]
    Ys2 <- Ys[jdx]
    modi <- fitSubtypeModel(Xs=Xs2, Ys=Ys2,
                            geneSet = geneSet,
                            breakVec=breakVec,
                            params=params,
                            caret.seed = caret.seed,
                            ptail=ptail)
    return(modi)
  }

  clusterExport(cl, 'Xs',  envir=environment())
  clusterExport(cl, 'Ys',  envir=environment())
  clusterExport(cl, 'geneSet',  envir=environment())
  clusterExport(cl, 'sampSize',  envir=environment())
  clusterExport(cl, 'seeds',  envir=environment())
  clusterExport(cl, 'caret.seed',  envir=environment())
  clusterExport(cl, 'breakVec',  envir=environment())
  clusterExport(cl, 'params',  envir=environment())
  clusterExport(cl, 'ptail',  envir=environment())
  clusterExport(cl, c('fitSubtypeModel','trainDataProc','cvFitOneModel','makeSetData','makeGenePairs','breakBin','binaryGene','featureSelection','testFun'),envir=environment())
  clusterExport(cl, c('xgb.DMatrix','xgb.cv','xgboost','trainControl','train','xgb.train'),envir=environment())
  clusterExport(cl, 'fitFun',  envir=environment())

  ens <- parLapply(cl=cl, X=1:n, fun = fitFun)

  res <- list(
    Repeat = list(
      Xs = Xs, Ys = Ys,
      geneSet = geneSet,
      n=n,
      sampSize=sampSize, sampSeed = sampSeed,
      breakVec=breakVec,
      params=params,
      ptail=ptail, numCores=numCores
    ),
    Model = ens
  )

  stopCluster(cl)

  return(res)
}








