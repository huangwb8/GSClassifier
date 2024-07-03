

#' @title geneMatch
#' @description Match the incoming data to what was used in training
#' @param X gene expression matrix, genes in rows, samples in columns
#' @param geneAnnotation A data frame with ENSEMBL, SYMBOL and ENTREZID of the genes in \code{geneSet}
#' @param geneid type of gene id in \code{X}. One of \code{'symbol'}, \code{'entrez'} and \code{'ensembl'}.
#' @param matchmode If your genes belong to ENSEMBL, SYMBOL or ENTREZ, you should set \code{matchmode = 'fix'}, in which all genes would be aligned to ENSEMBL. \code{'matchmode = free'} is a legacy feature, which means that you have to convert your genes into the same annotation type of the pre-trained \code{GSClassifier} model.
#' @import luckyBase
#' @return Matched index of genes to the expression matrix.
#' @export
geneMatch<- function(X,
                     geneAnnotation = NULL,
                     geneid='ensembl',
                     matchmode = c('fix', 'free')[1]){
  if(!matchmode %in% c('fix', 'free')){
    stop("Please input one of 'fix' and 'free'. ")
  } else if(matchmode == 'fix'){
    geneMatch_fixed(X,geneAnnotation,geneid)
  } else {
    geneMatch_free(X,geneAnnotation,geneid)
  }
}

#' @description fix mode of \code{geneMatch}
geneMatch_fixed <- function(X,
                      geneAnnotation = NULL,
                      geneid='ensembl') {

  ## Gene annotation
  if(is.null(geneAnnotation)){
    geneAnnotation <- readRDS(system.file("extdata", "general-gene-annotation.rds", package = "GSClassifier"))[["hg38"]] # Only for human
  }

  name_col <- c('SYMBOL','ENSEMBL','ENTREZID'); names(name_col) <- c('symbol','ensembl','entrez')

  # Single or Mutiple samples
  test <- is.matrix(X)|is.data.frame(X)

  if(test){

    # Multiple data
    if (!geneid %in% c('ensembl', 'symbol', 'entrez')){
      LuckyVerbose("For geneids, please use:  symbol, entrez, ensembl")
      return(NA)
    } else if (!geneid %in% 'ensembl'){
      rowname_ensembl <- convert(rownames(X), as.character(name_col[geneid]), 'ENSEMBL', geneAnnotation)
      X <- X[!is.na(rowname_ensembl),]
      rownames(X) <- convert(rownames(X), as.character(name_col[geneid]), 'ENSEMBL', geneAnnotation)
      # rowname_ensembl <- rowname_ensembl[!is.na(rowname_ensembl)]
      # X <- mergeMatrixDup(
      #   X,
      #   mergeCol = F,
      #   mergeRow = T,
      #   fun_row = mean,
      #   refRow = rowname_ensembl,
      #   verbose = F
      # )
      idx <- match(table = rownames(X), x = geneAnnotation$ENSEMBL)
      X2 <- X[idx,]
    } else {
      idx <- match(table = rownames(X), x = geneAnnotation$ENSEMBL)
      X2 <- X[idx,]
    }

      # Adds NA rows in missing genes
    rownames(X2) <- as.character(geneAnnotation[,name_col[names(name_col) == 'ensembl']])


  } else {

    # Single data

    if (!geneid %in% c('ensembl', 'symbol', 'entrez')){
      LuckyVerbose("For geneids, please use:  symbol, entrez, ensembl")
      return(NA)
    } else if (!geneid %in% 'ensembl'){
      name_ensembl <- convert(names(X), as.character(name_col[geneid]), 'ENSEMBL', geneAnnotation)
      X <- X[!is.na(name_ensembl)]
      names(X) <- convert(names(X), as.character(name_col[geneid]), 'ENSEMBL', geneAnnotation)
      idx <- match(table = names(X), x = geneAnnotation$ENSEMBL)
      X2 <- X[idx]
    } else {
      idx <- match(table = names(X), x = geneAnnotation$ENSEMBL)
      X2 <- X[idx]
    }

    # Adds NA rows in missing genes
    names(X2) <- as.character(geneAnnotation[,name_col[names(name_col) == 'ensembl']])

  }

  # Missing level
  matchError <- sum(is.na(idx)) / nrow(geneAnnotation)
  if(matchError >0){
    missGenes <- as.character(geneAnnotation[,name_col[names(name_col) == 'ensembl']])[is.na(idx)]
  } else {
    missGenes <- NA
  }

  # Output
  return(list(Subset=X2, matchError=matchError, missGenes = missGenes))
}

#' @description free mode of \code{geneMatch}
geneMatch_free<- function(X,
                      geneAnnotation = NULL,
                      geneid='ensembl') {

  ## Gene annotation
  if(is.null(geneAnnotation)){
    geneAnnotation <- readRDS(system.file("extdata", "modelGeneAnnotation.rds", package = "GSClassifier"))
  }

  name_col <- c('SYMBOL','ENSEMBL','ENTREZID'); names(name_col) <- c('symbol','ensembl','entrez')

  # Single or Mutiple samples
  test <- is.matrix(X)|is.data.frame(X)

  if(test){

    # Multiple data

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

    X2 <- X[idx,]  ### Adds NA rows in missing genes
    rownames(X2) <- as.character(geneAnnotation[,name_col[names(name_col) == geneid]])


  } else {

    # Single data

    if (geneid == 'symbol') {
      idx <- match(table = names(X), x = as.character(geneAnnotation$SYMBOL))  ### this is just for the EBPP genes ###

    } else if (geneid == 'entrez') {
      idx <- match(table = names(X), x = as.character(geneAnnotation$ENTREZID))

    } else if (geneid == 'ensembl') {
      idx <- match(table = names(X), x = geneAnnotation$ENSEMBL)

    } else {
      print("For geneids, please use:  symbol, entrez, ensembl")
      return(NA)
    }

    X2 <- X[idx]  ### Adds NA rows in missing genes
    names(X2) <- as.character(geneAnnotation[,name_col[names(name_col) == geneid]])

  }

  matchError <- sum(is.na(idx)) / nrow(geneAnnotation)

  return(list(Subset=X2, matchError=matchError))
}

#' @description Error report
#' @param geneMatchResult the \code{matchError} result of \code{\link{geneMatch}}
reportError <- function(geneMatchResult) {
  err <- geneMatchResult[['matchError']]
  LuckyVerbose(paste0("geneMatch: Percent missing genes=",err*100,"%."))
  if(err>0){
    LuckyVerbose('geneMatch: Missing genes=', paste0(geneMatchResult[['missGenes']], collapse = ', '))
  }
}

#' @description Correct X format
#' @inheritParams geneMatch
rightX <- function(X){

  if(F){
    X = expr2[,1]; names(X) <- rownames(expr2)
    X = as.matrix(expr2[,1]); rownames(X) <- rownames(expr2)
    # X = expr2[,1:2]
  }

  multi <- (is.matrix(X) | is.data.frame(X))

  if(multi){

    # Matrix

    if(ncol(X) == 1){
      # Single col matrix. To vector.
      rX <- as.vector(X); names(rX) <- rownames(X)
    } else {
      # Multiple matirx. Normal.
      rX <- X
    }

  } else {

    # Vector

    rX <- X

  }

  return(rX)

}



