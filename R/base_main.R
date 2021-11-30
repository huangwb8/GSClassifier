

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
#' @param err the \code{matchError} result of \code{\link{geneMatch}}
reportError <- function(err) {
  LuckyVerbose(paste0("Percent missing genes: ",err*100,"%."))
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



