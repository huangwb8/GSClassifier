
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
createPtail <-
  function(geneSet,
           interval = c(200, 400, 600, 800, 1000)) {
    # Test
    if (F) {
      ImmuneSubtype <-
        readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
      geneSet <- ImmuneSubtype$geneSet
      interval = c(200, 400, 600, 800, 1000)
    }

    Num_GeneSet <- length(geneSet)
    Num_Gene <- length(unique(unlist(geneSet)))

    ptails <- sapply(interval, function(x)solve_ptail(Num_GeneSet, Num_Gene, x))
    return(ptails/2)
  }

####=== Asistant fucntions ===####

solve_ptail <- function(Num_GeneSet, Num_Gene, interval) {

  # x = -b/(2a) ± √(b²-4ac) / (2a)

  a = Num_Gene ^ 2
  b = -(Num_Gene)
  c = - 2 * interval

  x1 <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
  # x2 <- (-b - sqrt(b^2 - 4*a*c))/(2*a)

  return(x1)
}
