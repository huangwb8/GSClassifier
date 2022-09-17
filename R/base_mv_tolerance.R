

#' @title mv_tolerance
#' @description Missing value tolerance analysis for \code{GSClassifier} models
#' @param X a data frame with the number of colume ≥ 2
#' @param gene.loss a vector of the number of genes masked
#' @param model a character of GSClassifier model or a \code{luckyModel} model
#' @param seed random seed
#' @param verbose whether to report the process
#' @importFrom pROC multiclass.roc
#' @return A list containing results of subtype identification and multi-roc analysis
#' @seealso \code{\link[pROC]{multiclass.roc}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' Empty
#' @export
mv_tolerance <- function(X,
                         gene.loss = c(2,4,6,8,10,12),
                         levels = c(1,2,3,4),
                         model = 'PAD.train_20220916',
                         seed = 487,
                         verbose = T){

  # Test
  if(F){

    # Internal validation cohort
    testData <- readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
    expr_pad <- testData$PanSTAD_expr_part
    modelInfo <- modelData(
      design = testData$PanSTAD_phenotype_part,
      id.col = "ID",
      variable = c("platform", "PAD_subtype"),
      Prop = 0.7,
      seed = 19871
    )
    validInform <- modelInfo$Data$Valid
    X <- expr_pad[,validInform$ID]

    # Other parameters
    gene.loss = c(2,4,6,8,10,12)
    model = 'PAD.train_20220916'
    seed = 487
    verbose = T
    levels = c(1,2,3,4)

  }

  # MVI with quantile algorithm
  X <- GSClassifier:::na_fill(X,
                              method='quantile',
                              seed = 411,
                              verbose = verbose)

  # Get model
  if(is.character(model)){
    # Model in GSClassifier
    m <- readRDS(system.file("extdata",
                             paste0(model,'.rds',collapse = ''),
                             package = "GSClassifier"))
  } else if(is.list(model)){
    # Model in luckyModel
    m <- model
  } else {
    # Error report
    stop('Please use a right model format.')
  }

  # Subtype identification for raw matrix
  res0 <- parCallEnsemble(
    X = X,
    ens = m$ens$Model,
    geneAnnotation = m$geneAnnotation,
    geneSet = m$geneSet,
    scaller = m$scaller$Model,
    geneid = "ensembl",
    subtype = NULL,
    verbose = verbose,
    numCores = 4
  )

  # multi-ROC: time-consuming if a large matrix used
  set.seed(seed); seeds <- sample(1:100000,
                                  length(gene.loss),
                                  replace = F)
  mAUC <- list(); model_res <- list()
  for(i in 1:length(gene.loss)){ # i=1

    # Masked matrix with zero value
    X2 <- masked_matrix(X,
                        gene.loss = gene.loss[i],
                        seed = 487)

    # Subtype identification
    res <- parCallEnsemble(
      X = X2,
      ens = m$ens$Model,
      geneAnnotation = m$geneAnnotation,
      geneSet = m$geneSet,
      scaller = m$scaller$Model,
      geneid = "ensembl",
      subtype = NULL,
      verbose = verbose,
      numCores = 4
    )

    mAUC[[paste0('GeneLoss=', gene.loss[i])]] <- multiclass.roc(
      response = res0$BestCall,
      predictor = res$BestCall,
      levels = levels,
      quiet = T)
    model_res[[paste0('GeneLoss=', gene.loss[i])]] <- res
  }
  # mymusic()

  # Output
  l <- list(
    multiAUC = mAUC,
    modelRes = model_res
  )
  return(l)

  # saveRDS(l, './data/results-mv_tolerance-GC.rds')

}

#' @inheritParams mv_tolerance
masked_matrix <- function(X,
                          gene.loss = 2,
                          seed = 487){

  # Zero Matrix
  # Xzero <- matrix(
  #   data = rep(0, nrow(X)*ncol(X)),
  #   nrow = nrow(X), byrow = F,
  #   dimnames = list(rownames(X),
  #                   colnames(X)))

  # Random seeds
  set.seed(seed); seeds <- sample(1:100000, ncol(X), replace = F)
  X2 <- as.matrix(X)
  for(i in 1:ncol(X)){ # i=1
    set.seed(seeds[i])
    X2[sample(1:nrow(X2), gene.loss, replace = F),i] <- 0
  }
  return(X2)
}



