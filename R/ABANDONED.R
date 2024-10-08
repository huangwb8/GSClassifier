
# Legacy codes

if (F) {
  #' @description Compute pairwise set comparisons for gene expression data.
  #' @param Xmat Matrix of gene expression data.
  #' @param geneSet List of genes for classification.
  #' @return Matrix of pairwise set comparison results.
  #' @details
  #' A version by Claude-3.5-sonnet
  makeSetData_legacy_02 <- function(Xmat, geneSet) {
    nGS <- length(geneSet)

    # Pre-compute feature names
    featureNames <-
      combn(
        nGS,
        2,
        FUN = function(x)
          paste0('s', x[1], 's', x[2])
      )

    # Pre-compute gene set indices
    geneSetIndices <-
      lapply(geneSet, function(set)
        which(rownames(Xmat) %in% set))

    # Compute results for all samples at once
    resultMat <- vapply(1:ncol(Xmat), function(i) {
      # i=1
      sapply(1:length(featureNames), function(idx) {
        # idx=1
        setNum <- Fastextra(featureNames[idx], 's') %>% as.numeric()
        j1 <- setNum[2]
        j2 <- setNum[3]

        vals1 <- Xmat[geneSetIndices[[j1]], i]
        vals2 <- Xmat[geneSetIndices[[j2]], i]

        if (length(vals1) == 0 | length(vals2) == 0) {
          return(0)
        } else {
          comparisons <- outer(vals1, vals2, ">")
          return(sum(comparisons) / length(comparisons))
        }
      })
    }, numeric(choose(nGS, 2)))

    rownames(resultMat) <- featureNames
    colnames(resultMat) <- colnames(Xmat)

    return(resultMat)
  }

  # Bensz version
  makeSetData_legacy_01 <- function(Xmat, geneSet) {
    resultList <- list()
    nGS <- length(geneSet) # nGS=41

    featureNames <- c()
    for (j1 in 1:(nGS - 1)) {
      for (j2 in (j1 + 1):nGS) {
        featureNames <- c(featureNames, paste0('s', j1, 's', j2))
      }
    }

    # for each sample
    for (i in 1:ncol(Xmat)) {
      res0 <- numeric(length = length(featureNames))
      idx <- 1
      for (j1 in 1:(nGS - 1)) {
        for (j2 in (j1 + 1):nGS) {
          set1 <- geneSet[[j1]]
          set2 <- geneSet[[j2]]
          vals1 <- Xmat[rownames(Xmat) %in% set1, i]
          vals2 <- Xmat[rownames(Xmat) %in% set2, i]
          if (length(vals1) == 0 | length(vals2) == 0) {
            res0[idx] <- 0 # weight it as 0
          } else {
            res1 <- sapply(vals1, function(v1)
              sum(v1 > vals2, na.rm = T))
            res0[idx] <-
              sum(res1, na.rm = T) / (length(vals1) * length(vals2))
          }
          idx <- idx + 1
        }
      }
      resultList[[i]] <- as.numeric(res0)
    }
    resMat <- do.call(cbind, resultList)
    colnames(resMat) <- colnames(Xmat)
    rownames(resMat) <- featureNames
    return(resMat)
  }

  makeGenePairs_legacy_before20241008 <- function(genes, Xsub) {
    # collect results here
    resList <- list()
    # all pairs of genes here
    all.combos <- t(combn(genes, 2))

    for (i in 1:nrow(all.combos)) {
      gi <- all.combos[i, 1]  # the pivot gene
      gval <- as.numeric(Xsub[gi, ]) # and it's value
      # get the paired genes out.
      gcom <- all.combos[all.combos[, 1] == gi, 2]
      # pick sample j, get pivot value j which is for gene i, and pivot all the paired genes
      res0 <-
        lapply(1:ncol(Xsub), function(j)
          binaryGene(gval[j], Xsub[gcom, j]))
      # make matrix of features.
      resMat <- do.call('cbind', res0)
      rownames(resMat) <- sapply(gcom, function(gj)
        paste0(gi, ':', gj))
      resList[[gi]] <- resMat
    }
    Xbin <- do.call('rbind', resList)
    colnames(Xbin) <- colnames(Xsub)
    return(Xbin)
  }

  # cvFitOneModel
  cvFitOneModel_legacy <- function(Xbin,
                                   Ybin,
                                   params = list(
                                     # xgb.cv only
                                     nfold = 5,
                                     nrounds = 100,
                                     # xgb.cv & xgboost
                                     max_depth = 10,
                                     eta = 0.5,
                                     nthread = 5,
                                     colsample_bytree = 1,
                                     min_child_weight = 1
                                   ),
                                   breakVec = c(0, 0.25, 0.5, 0.75, 1.0),
                                   genes,
                                   verbose = F) {
    # Test
    if (F) {
      # Example
      library(GSClassifier)
      library(xgboost)
      testData <-
        readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
      expr <- testData$PanSTAD_expr_part
      design <- testData$PanSTAD_phenotype_part
      modelInfo <- modelData(
        design,
        id.col = "ID",
        variable = c("platform", "PAD_subtype"),
        Prop = 0.1,
        seed = 145
      )
      Xs <- expr[, modelInfo$Data$Train$ID]
      y <- modelInfo$Data$Train
      y <- y[colnames(Xs), ]
      Ys <-
        ifelse(y$PAD_subtype == 'PAD-I',
               1,
               ifelse(
                 y$PAD_subtype == 'PAD-II',
                 2,
                 ifelse(
                   y$PAD_subtype == 'PAD-III',
                   3,
                   ifelse(y$PAD_subtype == 'PAD-IV', 4, NA)
                 )
               ))
      table(Ys) / length(Ys)
      PADi <-
        readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
      geneSet <- PADi$geneSet
      res <- trainDataProc(
        Xmat = Xs,
        Yvec = Ys,
        geneSet = geneSet,
        subtype = 2,
        ptail = 0.5,
        breakVec = c(0, 0.25, 0.5, 0.75, 1)
      )

      # Data
      Xbin <- res$dat$Xbin
      Ybin <- res$dat$Ybin
      genes <- res$dat$Genes
      params <- list(
        # xgboost & xgb.cv
        nfold = 5,
        nrounds = 100,

        # xgboost
        max_depth = 10,
        colsample_bytree = 1,
        min_child_weight = 1,
        eta = 0.5,
        gamma = 0.25,
        subsample = 0.7
      )
      breakVec = c(0, 0.25, 0.5, 0.75, 1.0)
      verbose = T

    }

    dtrain <- xgb.DMatrix(Xbin, label = Ybin)

    # xgb.cv

    # 2022-09-15  : WARNING: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.

    # 2022-09-16: WARNING: If you are loading a serialized model (like pickle in Python, RDS in R) generated by
    # older XGBoost, please export the model by calling `Booster.save_model` from that version
    # first, then load it back in current version. See: https://xgboost.readthedocs.io/en/latest/tutorials/saving_model.html
    # for more details about differences between saving model and serializing.

    for (i in 1:10000) {
      x <- tryCatch(
        cvRes <- xgb.cv(
          params = params,
          nrounds = params$nrounds,
          nfold = params$nfold,
          data = dtrain,
          early_stopping_rounds = 2,
          metrics = list("logloss", "auc"),
          objective = "binary:logistic",
          verbose = verbose
        ),
        error = function(e)
          e)
      if ('message' %in% names(x)) {
        if (verbose)
          LuckyVerbose('Attention! AUC: the dataset only contains pos or neg samples. Repeat xgb.cv')
        x_error <- x
      } else {
        cvRes <- x
        break
      }
    }

    if (verbose)
      LuckyVerbose('Best interation: ', cvRes$best_iteration)

    # xgboost via best interation
    bst <- xgboost(
      params = params,
      data = Xbin,
      label = Ybin,
      nrounds = cvRes$best_iteration,
      objective = "binary:logistic",
      verbose = ifelse(verbose, 1, 0)
    )

    return(list(
      bst = bst,
      breakVec = breakVec,
      genes = genes
    ))
  }



  # nround不是xgb.cv选的。容易过拟合。
  cvFitOneModel <- function(Xbin,
                            Ybin,
                            genes,
                            params = list(
                              # xgb.cv only
                              nfold = 5,
                              nrounds = 15,
                              # xgb.cv & xgboost
                              max_depth = 10,
                              eta = 0.5,
                              nthread = 5,
                              colsample_bytree = 1,
                              min_child_weight = 1
                            ),
                            breakVec = c(0, 0.25, 0.5, 0.75, 1.0),
                            seed = 102,
                            verbose = F) {
    # Test
    if (F) {
      # Example
      library(GSClassifier)
      library(xgboost)
      testData <-
        readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
      expr <- testData$PanSTAD_expr_part
      design <- testData$PanSTAD_phenotype_part
      modelInfo <- modelData(
        design,
        id.col = "ID",
        variable = c("platform", "PAD_subtype"),
        Prop = 0.1,
        seed = 145
      )
      Xs <- expr[, modelInfo$Data$Train$ID]
      y <- modelInfo$Data$Train
      y <- y[colnames(Xs), ]
      Ys <-
        ifelse(y$PAD_subtype == 'PAD-I',
               1,
               ifelse(
                 y$PAD_subtype == 'PAD-II',
                 2,
                 ifelse(
                   y$PAD_subtype == 'PAD-III',
                   3,
                   ifelse(y$PAD_subtype == 'PAD-IV', 4, NA)
                 )
               ))
      table(Ys) / length(Ys)
      PADi <-
        readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
      geneSet <- PADi$geneSet
      res <- trainDataProc(
        Xmat = Xs,
        Yvec = Ys,
        geneSet = geneSet,
        subtype = 2,
        ptail = 0.5,
        breakVec = c(0, 0.25, 0.5, 0.75, 1)
      )

      # Data
      Xbin <- res$dat$Xbin
      Ybin <- res$dat$Ybin
      genes <- res$dat$Genes
      params <- list(
        # xgboost & xgb.cv
        nrounds = 15,

        # xgboost
        max_depth = 10,
        colsample_bytree = 1,
        min_child_weight = 1,
        eta = 0.5,
        gamma = 0.25,
        subsample = 0.7
      )
      breakVec = c(0, 0.25, 0.5, 0.75, 1.0)
      verbose = T
      seed = 102

    }

    # xgboost via best interation
    set.seed(seed)
    params_xg <- params[-match(c('nrounds'), names(params))]
    bst <- xgboost(
      params = params_xg,
      data = Xbin,
      label = Ybin,
      objective = "binary:logistic",
      nrounds = params$nrounds,
      verbose = ifelse(verbose, 1, 0)
    )
    return(list(
      bst = bst,
      breakVec = breakVec,
      genes = genes
    ))
  }

}
