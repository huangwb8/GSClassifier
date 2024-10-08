
####====================Basic function==================####

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
  for (gi in genePairs) {
    pairList[[gi[1]]] <- c(pairList[[gi[1]]], gi[2])
  }

  resList <- list() # then for each pivot gene
  for (gi in names(pairList)) {
    # assuming it's in the data ... really should be!
    if (gi %in% rownames(X)) {
      gs <-
        pairList[[gi]]                  ## can end up with the pivot gene in the genes...
      pval <-
        as.numeric(X[gi, ])            ## pivot values across samples
      idx <-
        match(table = rownames(X), x = gs) ## get index to genes for this pivot
      Xsub <-
        X[idx, ]                       ## subset the matrix, NAs for missing genes, pivot gene on top
      is_num_vector <-
        class(Xsub)[1]  %in% c('numeric', "integer") & length(gs) == 1
      if (is_num_vector) {
        Xsub <- matrix(data = Xsub,
                       ncol = ncol(X),
                       nrow = 1)
        colnames(Xsub) <- colnames(X)
      }
      rownames(Xsub) <- gs                  ## give gene IDs
      res0 <-
        lapply(1:ncol(Xsub), function(a)
          binaryGene(pval[a], Xsub[, a]))  ## create binary values
      resList[[gi]] <- do.call('cbind', res0)
    } else {
      # else we need to include some dummy rows
      randMat <-
        matrix(data = rbinom(n = length(pairList[[gi]]) * ncol(X), prob = 0.5, 1),
               ncol = ncol(X))
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
spliteMatrix <- function(X, cutoff = 20) {
  # Cut_vector
  cut_vector <- function(vt, nsplit = 100) {
    if (nsplit == 1) {
      v2 <- list(vt)

      names(v2) <- paste0("1-", length(v2))
    } else {
      ##转化成位置向量
      len.vt = 1:length(vt)

      ##间隔
      len1 <- floor(length(len.vt) / nsplit)

      ## low.ci and upper.ci
      low.ci <- 1

      for (i in 1:(nsplit - 1)) {
        low.ci[i + 1] <- low.ci[i] + len1
      }
      upper.ci <- len1
      for (i in 1:(nsplit - 1)) {
        upper.ci[i + 1] <- upper.ci[i] + len1
      }

      ## upper.ci的最后一个值
      upper.ci[length(upper.ci)] <- length(vt)

      ## 形成列表
      v2 <- NULL
      for (i in 1:length(low.ci)) {
        v2.i <- low.ci[i]:upper.ci[i]
        v2.i <- vt[v2.i]
        v2 <- c(v2, list(v2.i))
        names(v2)[i] <- paste0(low.ci[i], "-", upper.ci[i])
      }
    }
    ## 输出结果
    return(v2)
  }

  # Splite matrix
  nXL <- floor(ncol(X) / cutoff)
  eXL <- cut_vector(colnames(X), nXL)
  XL <- list()
  for (i in 1:length(eXL)) {
    # i=1
    XL[[names(eXL)[i]]] <- X[, eXL[[i]]]
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
                     mods = NULL,
                     geneSet,
                     nClust = 4) {
  # working with matrices
  Xmat <- as.matrix(X)

  # NA filling with recursive partitioning and regression trees
  # Xmat <- na_fill(Xmat, method="anova", na.action = na.rpart)

  # Set features
  nGS <- length(geneSet)
  featureNames <- c()
  for (j1 in 1:(nGS - 1)) {
    for (j2 in (j1 + 1):nGS) {
      featureNames <- c(featureNames, paste0('s', j1, 's', j2))
    }
  }

  # get out the relevant items
  breakVec <- mods$breakVec
  genes <- mods$bst$feature_names
  singleGenes <- genes[!str_detect(genes, ':')]
  singleGenes <- singleGenes[!singleGenes %in% featureNames]
  pairedGenes <- genes[str_detect(genes, ':')]
  setFeatures <- genes[str_detect(genes, 's[0-9]{1,6}s[0-9]{1,6}')]

  # bin the expression data
  Xbinned <- apply(Xmat, 2, breakBin, breakVec)
  rownames(Xbinned) <- rownames(Xmat)

  # and subset the genes to those not in pairs
  if (length(singleGenes) > 1) {
    Xbinned <- Xbinned[singleGenes, ]
  }

  # here we have expression data, and we're using the pairs model
  # so we need to make pairs features.
  Xpairs <- createPairsFeatures(Xmat, pairedGenes)
  colnames(Xpairs) <- colnames(Xmat)

  # gene set features.
  Xset <- makeSetData(Xmat, geneSet)
  if (length(setFeatures) > 1) {
    Xset <- Xset[setFeatures, ]
  }

  # join the data types and transpose
  Xbin <- t(rbind(Xbinned, Xpairs, Xset))
  return(Xbin)
}


#' @description Make subtype calls for one sample
#' @param mods xgboost model list
#' @param X gene expression matrix, genes in rows, samples in columns
#' @param ci cluster label, and index into mods
#' @param verbose whether report messages
#' @inheritParams dataProc
#' @import xgboost
#' @return preds of one cluster model.
callOneSubtype <-
  function(mods, X, ci, geneSet, nClust, verbose = T) {
    # 2023-12-26：List in R is special. If the name of the list is 1,2,3,4, when you call l[[1]], it would call the name of "1" instead of the first one.
    # Xbin needs to have the same columns as the training matrix.

    # if(verbose) (paste0('calling subtype ', ci))

    if (as.character(ci) %in% names(mods)) {
      if (length(mods[[as.character(ci)]]) > 1) {
        mi <- mods[[as.character(ci)]] # or we can use: mi <- mods[[ci]]
        Xbin <- dataProc(X, mods = mi, geneSet, nClust)

        # WARNING: amalgamation/../src/learner.cc:438:
        # If you are loading a serialized model (like pickle in Python, RDS in R) generated by
        # older XGBoost, please export the model by calling `Booster.save_model` from that version
        # first, then load it back in current version. See:
        #   https://xgboost.readthedocs.io/en/latest/tutorials/saving_model.html
        # for more details about differences between saving model and serializing.

        # Prediction
        pred <- predict(mi$bst, Xbin)
      } else {
        pred <- rep(0, ncol(X))
      }
    } else {
      # if(verbose) paste0('calling subtype ', ci, ': No avaliable model. Return 0!')

      # The model is unavailable. Return 0.
      pred <- rep(0, ncol(X))

    }

    return(pred)
  }


#' @description Make subtype calls for each sample
#' @inheritParams callOneSubtype
callSubtypes <- function(mods, X, geneSet, clusterName, verbose = T) {
  pList <-
    lapply(clusterName, function(mi)
      callOneSubtype(mods, X, mi, geneSet, nClust, verbose))
  pMat  <- do.call('cbind', pList)
  colnames(pMat) <- clusterName # names(mods)
  bestCall <-
    apply(pMat, 1, function(pi)
      colnames(pMat)[which(pi == max(pi))][1])
  # 2022-9-27. There's no problem because bestcall_max is a legacy category strategy and would not be recommended.
  # bestCall <- apply(pMat, 1, function(pi) colnames(pMat)[which(pi == max(pi)[1])])
  #     1           2           3           4
  # 0.077546641 0.004749997 0.832491577 0.832491577
  return(data.frame(
    SampleID = colnames(X),
    BestCall = bestCall,
    pMat,
    stringsAsFactors = F
  ))
}
