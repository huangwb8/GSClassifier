

#' @title PADi
#' @description PAD for individuals
#' @param cancer.type The type of cancer
#' @param version model version
#' @inheritParams parCallEnsemble
#' @return A data frame containing personalized PAD subtypes
#' @seealso \code{\link{callEnsemble}}; \code{\link{parCallEnsemble}}.
#' @details The current software only supports parameters: cancer.type = "GC", version = "20200110". More cancer types or model version were under developing.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' # Package
#' if (!requireNamespace("GSClassifier", quietly = TRUE)){
#'   devtools::install_github("huangwb8/GSClassifier")
#' }
#' library(GSClassifier)
#'
#' # Load test data of RNA-Seq expression
#' testData <- readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))[['Kim2018_3']]
#' dim(testData) # 19118 3
#'
#' # Data formats (Different data but the same function and usage)
#'
#' # 1.Mutiple data
#' X = testData
#'
#' # 2.Single data
#' # X <- testData[,1]; names(X) <- rownames(testData)
#' # X <- as.matrix(testData[,1]); rownames(X) <- rownames(testData)
#'
#' ## Method1: use a general function called 'callEnsemble'
#' res_padi <- callEnsemble(
#'   X = X,
#'   ens = NULL,
#'   geneAnnotation = NULL,
#'   geneSet = NULL,
#'   scaller = NULL,
#'   geneid = "ensembl",
#'   matchmode = 'fix',
#'   subtype = "PAD.train_20200110",
#'   verbose = F
#' )
#'
#' ## Method2: use a specific function called 'PADi'
#' res_padi <- PADi(X = X, verbose = F)
#'
#' ## Method3: Parallel strategy for lots of samples (empirically >50) to save time. Not Run for small cohorts.
#' # res_padi <- parCallEnsemble(
#' #   X = X,
#' #   ens = NULL,
#' #   geneAnnotation = NULL,
#' #   geneSet = NULL,
#' #   scaller = NULL,
#' #   geneid = 'ensembl',
#' #   matchmode = 'fix',
#' #   subtype = 'PAD.train_20200110',
#' #   verbose = T,
#' #   numCores = 2)
#' # res_padi <- PADi(X = X, verbose = F, numCores = 4)
#' @export
PADi <- function(X,
                 geneid = "ensembl",
                 matchmode = c('fix', 'free')[1],
                 cancer.type = 'GC',
                 version = c('20200110',
                             '20220916')[1],
                 numCores = 0,
                 verbose = T) {
  # Test
  if (F) {
    X = expr2
    geneid = "ensembl"
  }


  pass <- (cancer.type == 'GC')

  subtype <- paste0("PAD.train_", version)

  if (pass & numCores == 0) {
    # Common mode
    callEnsemble(
      X = X,
      ens = NULL,
      geneAnnotation = NULL,
      geneSet = NULL,
      scaller = NULL,
      geneid = geneid,
      matchmode = matchmode,
      # subtype = "PAD.train_20200110",
      subtype = subtype,
      verbose = verbose
    )

  } else if (pass & numCores > 0) {
    # Parallel mode
    parCallEnsemble(
      X = X,
      ens = NULL,
      geneAnnotation = NULL,
      geneSet = NULL,
      scaller = NULL,
      geneid = geneid,
      matchmode = matchmode,
      subtype = subtype,
      verbose = verbose,
      numCores = numCores
    )

  } else{
    LuckyVerbose(
      "Incorrect 'cancer.type' or 'version'. Only cancer.type = 'GC' is available in the current software."
    )

  }

}
