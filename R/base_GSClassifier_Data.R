


#' @title Data in \code{GSClassifier} package
#' @description Internal dataset of GSClassifier package
#' @param ImmuneSubtype Data integration for calling of immune subtype by Thorsson et al. in 2018.
#' @param PAD Data integration for calling of pan-immune activation/dysfunction subtype (PAD).
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' PAD <- readRDS(system.file("extdata", "PAD.train_20200110", package = "GSClassifier"))
#' ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
#' @export
GSClassifier_Data <- function(ImmuneSubtype=NULL,
                              PAD=NULL){

  data.p <- system.file("extdata", package = "GSClassifier")
  data.n <- list.files(data.p,pattern = '.rds$',full.names = F,recursive = F)

  ## Available data
  message('Available data:')
  for(i in data.n) cat(' ',i,'\n')

  ## Example
  message('Usage example: ')
  cat('  PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier"))','\n')
  cat('  ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))')

}


