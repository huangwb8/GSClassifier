
#' @title Greate Ptail Vector
#' @description Create Ptail Vector
#' @param interval Numeric. Define the number of features you want in model training.
#' @inheritParams subtypeVector
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
#' geneSet <- ImmuneSubtype$geneSet
#' ptail <- createPtail(geneSet, interval = c(200,400,600,800,1000))
#' print(ptail)
#' @export
createPtail <- function(geneSet, interval = c(200,400,600,800,1000)){

  # Test
  if(F){
    ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
    geneSet <- ImmuneSubtype$geneSet
    interval = c(200,400,600,800,1000)
  }

  Num_GeneSet <- length(geneSet)
  Num_Gene <- length(unique(unlist(geneSet)))
  Num_Feature <- Num_GeneSet * (Num_GeneSet - 1)/2 + Num_Gene * (Num_Gene - 1)/2

  interval <- interval[interval<=Num_Feature]
  if(length(interval) == 0){
    stop('Please use lower interval than ', Num_Feature, '!')
  } else {
    return(interval/Num_Feature)
  }
}
