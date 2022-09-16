


####====================Basic function==================####

#' @title NA filling with recursive partitioning and regression trees
#' @description NA filling with recursive partitioning and regression trees
#' @param Xmat A gene expression matrix with NA value
#' @param method One of \code{'quantile'}, \code{'rpart'} and \code{NULL}.
#' @importFrom luckyBase is.one.na
#' @import rpart
#' @seealso \code{\link[rpart]{rpart}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
na_fill <- function(Xmat,
                    method=c('quantile','rpart',NULL)[1],
                    seed=2022,
                    verbose=T){

  if(is.null(method)){

    if(verbose) LuckyVerbose('Without missing value imputation... ')
    Xmat

  } else if(method == 'rpart'){
    if(verbose) LuckyVerbose('Missing value imputation with RPART algorithm...')
    na_fill_rpart(Xmat,
                  method="anova",
                  na.action = na.rpart)

  } else if(method == 'quantile'){
    if(verbose) LuckyVerbose('Missing value imputation with quantile algorithm...')
    na_fill_quantile(Xmat,seed)

  } else {
    if(verbose) LuckyVerbose('Without missing value imputation... ')
    Xmat
  }
}

#' @inheritParams rpart::rpart
#' @inheritParams na_fill
na_fill_rpart <- function(Xmat,
                          method="anova",
                          na.action = na.omit){

  Xmat_2 <- t(Xmat)
  na.pos <- apply(Xmat_2,2,is.one.na)

  ## Test whether there're some NAs
  if(T %in% na.pos){

    # With NA value
    Xmat_2 <- as.data.frame(Xmat_2)
    na.marker <- names(na.pos)[na.pos]
    for(i in 1:length(na.marker)){ # i=1
      m.i <- na.marker[i] # m.i
      f.i <- as.formula(paste0(m.i," ~ ."))
      Xmat.i <- Xmat_2[!is.na(Xmat_2[,m.i]),]
      Xmat.i.na <- Xmat_2[is.na(Xmat_2[,m.i]),]
      anova_mod <- rpart(f.i, data=Xmat.i, method=method, na.action=na.action)
      anova_pred <- predict(anova_mod,Xmat.i.na) # View(Xmat_2[is.na(Xmat_2[,m.i]),])
      Xmat_2[rownames(Xmat.i.na),m.i] <- anova_pred
    }
  }

  # anyNA(Xmat_2) # FALSE
  return(t(Xmat_2))
  # ?rpart::rpart
  # ?base::anyNA
}

#' @inheritParams na_fill
na_fill_quantile <- function(Xmat,
                             seed=2022){

  # Test
  if(F){
    testData <- readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
    Xmat <- testData$PanSTAD_expr_part
  }

  na.pos <- apply(Xmat,2,is.one.na)
  set.seed(seed); seeds <- sample(1:ncol(Xmat)*10, sum(na.pos), replace = F)
  tSample <- names(na.pos)[na.pos]
  quantile_vector <- (1:1000)/1000

  for(i in 1:length(tSample)){ # i=1

    sample.i <- tSample[i]
    expr.i <- Xmat[, sample.i]
    expr.i.max <- max(expr.i, na.rm = T)
    expr.i.min <- min(expr.i, na.rm = T)
    set.seed(seeds[i]);
    expr.i[is.na(expr.i)] <-
      expr.i.min +
      (expr.i.max-expr.i.min) * sample(quantile_vector,
                                       sum(is.na(expr.i)),
                                       replace = T)
    Xmat[, sample.i] <- expr.i
  }

  return(Xmat)

}


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
                             ptail = 0.5) {
  idx <- which((testRes < quantile(testRes, ptail, na.rm = T)) |
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
  nGS <- length(geneSet) # nGS=4

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

  # NA filling with recursive partitioning and regression trees
  # Attention! This would make the model training more time-consuming. This step could be done in the begining of GSClassfier::fitSubtypeModel
  # Xmat <- na_fill(Xmat, method="anova", na.action = na.rpart)

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec) # bin each column
  rownames(Xbinned) <- rownames(Xmat)

  # rank the data, and use it for feature selection
  Xrank <- apply(Xmat, 2, rank)
  testRes <- sapply(1:nrow(Xrank), function(gi) testFun(as.numeric(Xrank[gi,]), Ybin))  # get genes with big rank diffs.

  # subset the expression data for pairs
  Xfeat <- featureSelection(Xmat, Ybin,
                            testRes=testRes,
                            ptail=ptail)  # subset genes
  Xpairs <- makeGenePairs(Xfeat$Genes, Xfeat$Xsub)

  # subset the binned genes
  Xbinned <- Xbinned[Xfeat$Genes,]

  # gene set features.
  nGS <- length(geneSet) # Optimized for single signature
  if(nGS == 1){
    # Without gene sets interaction
    Xset <- NULL
  } else {
    # With gene sets interaction
    Xset <- makeSetData(Xmat,geneSet) #--bug--
  }

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
#' @param verbose whether report modeling process
#' @inheritParams trainDataProc
#' @importFrom xgboost xgboost xgb.cv xgb.DMatrix
#' @return A single xgboost classifier.
cvFitOneModel <- function(Xbin, Ybin,
                          params=list(max_depth = 2,
                                      eta = 0.5,
                                      nrounds = 100,
                                      nthread = 5, nfold=5),
                          breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                          genes,
                          verbose = F){
  dtrain <- xgb.DMatrix(Xbin, label = Ybin)

  # xgb.cv

  # 2022-09-15  : WARNING: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.

  # 2022-09-16: WARNING: If you are loading a serialized model (like pickle in Python, RDS in R) generated by
  # older XGBoost, please export the model by calling `Booster.save_model` from that version
  # first, then load it back in current version. See: https://xgboost.readthedocs.io/en/latest/tutorials/saving_model.html
  # for more details about differences between saving model and serializing.

  for(i in 1:10000){
    x <- tryCatch(
      cvRes <- xgb.cv(data = dtrain,
                      nrounds=params$nrounds,
                      nthread=params$nthread,
                      nfold=params$nfold,
                      max_depth=params$max_depth,
                      eta=params$eta,
                      early_stopping_rounds=2,
                      metrics = list("logloss", "auc"), # 'logloss' instead of 'error'
                      objective = "binary:logistic",
                      verbose = verbose),
      error = function(e)e)
    if('message' %in% names(x)){
      if(verbose) LuckyVerbose('Attention! AUC: the dataset only contains pos or neg samples. Repeat xgb.cv')
      x_error <- x
    } else {
      cvRes <- x
      break
    }
  }

  if(verbose) LuckyVerbose('Best interation: ',cvRes$best_iteration)

  # xgboost via best interation
  bst <- xgboost(data = Xbin,
                 label = Ybin,
                 max_depth=params$max_depth,
                 eta=params$eta,
                 nrounds = cvRes$best_iteration,
                 nthread=params$nthread,
                 objective = "binary:logistic",
                 verbose = ifelse(verbose,1,0))

  return(list(bst=bst, breakVec=breakVec, genes=genes))
}

#' @description Train a single subtype model using cross validation
#' @importFrom xgboost xgb.train xgb.DMatrix
#' @importFrom caret trainControl train
#' @inheritParams cvFitOneModel
cvFitOneModel2 <- function(Xbin, Ybin,
                           breakVec=c(0, 0.25, 0.5, 0.75, 1.0),
                           genes,
                           seed = 101,
                           caret.grid = NULL,
                           verbose = F){

  # Data
  x = Xbin
  y = factor(Ybin,levels = c(0,1))

  # Parameters of caret::train
  grid <- caret.grid

  cntrl <- trainControl(
    method = "cv",
    number = 5,
    verboseIter = verbose,
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
                       nrounds = param$nrounds,
                       verbose = ifelse(verbose,1,0))

  return(list(bst=xgb.fit, breakVec=breakVec, genes=genes))


}


