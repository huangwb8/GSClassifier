
#' @title Data in \code{GSClassifier} package
#' @description Internal dataset of GSClassifier package
#' @param model logic. Whether to only show a list of \code{GSClassifier} models
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @details ImmuneSubtype is a model for personalized calling of immune subtype by Thorsson et al. in 2018. \cr
#' PAD.train_20200110 is a model for personalized calling of Pan-immune Activation Dysfunction Subtypes(PAD).
#' @examples
#' PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier"))
#' ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
#' @export
GSClassifier_Data <- function(model = T) {
  data.p <- system.file("extdata", package = "GSClassifier")
  data.n <-
    list.files(
      data.p,
      pattern = '.rds$',
      full.names = F,
      recursive = F
    )
  if (model)
    data.n <-
    data.n[!data.n %in% c('general-gene-annotation.rds', 'testData.rds')]

  ## Available data
  message('Available data:')
  for (i in data.n)
    cat(' ', i, '\n')

  ## Example
  message('Usage example: ')
  cat(
    '  PAD <- readRDS(system.file("extdata", "PAD.train_20200110.rds", package = "GSClassifier"))',
    '\n'
  )
  cat(
    '  ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))'
  )

}
