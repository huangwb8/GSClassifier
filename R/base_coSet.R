
#' @rdname CCS-function.coSet
#' @title coSet
#' @description Estimate co-expression genes distribution
#' @param sets a list containing series of gene symbols
#' @inheritParams ccs
#' @importFrom luckyBase LuckyVerbose
#' @import foreach
#' @import parallel
#' @import doParallel
#' @return a count matrix
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' set.seed(4090)
#' elements <- 1:6
#' sets <- list(
#'   sample(elements, 3),
#'   sample(elements, 4),
#'   sample(elements, 2),
#'   sample(elements, 5),
#'   sample(elements, 3))
#' resM <- coSet(sets)
#' print(resM)
#' @export
coSet <- function(sets,
                  verbose = T,
                  numCores = 1) {
  # Test
  if (F) {
    library(foreach)
    library(parallel)
    library(doParallel)
    sets <- list(
      sample(elements, 3),
      sample(elements, 4),
      sample(elements, 2),
      sample(elements, 5),
      sample(elements, 3)
    )
  }

  # Element
  element <- unique(unlist(sets))
  c <- length(element)

  # Empty matrix
  M <-
    matrix(
      0,
      nrow = c,
      ncol = c,
      dimnames = list(element, element)
    )

  # Add value
  time_int_1 <- system.time({
    if (numCores <= 1) {
      for (i in 1:length(sets)) {
        # i=1
        if (verbose)
          LuckyVerbose('coSet: processing GeneSet - ', names(sets)[i], '...')
        set <- sets[[i]]
        c_set <- length(set)
        M1 <-
          matrix(
            1,
            nrow = c_set,
            ncol = c_set,
            dimnames = list(set, set)
          )
        M <- M + reshapeMatrix(M1, M)
      }
    } else {
      # Use parallel strategy
      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      M <- foreach(i = 1:length(sets), .combine = '+') %dopar% {
        set <- sets[[i]]
        c_set <- length(set)
        M1 <-
          matrix(
            1,
            nrow = c_set,
            ncol = c_set,
            dimnames = list(set, set)
          )
        GSClassifier:::reshapeMatrix(M1, M)
      }
      stopCluster(cl)
    }
  })
  if (verbose)
    LuckyVerbose('coSet: Comsuming time = ', round(sum(time_int_1, na.rm = T), 2), 's...')

  # Output
  if (verbose)
    LuckyVerbose('coSet: All done!')
  return(M)

}


#### Assistant functions ####
reshapeMatrix <- function(minM, maxM) {
  # Test
  if (F) {
    minM = M1
    maxM = M
  }

  # Check
  check <-
    identical(rownames(minM), colnames(minM)) &
    identical(rownames(maxM), colnames(maxM)) &
    all(rownames(minM) %in% rownames(maxM))
  if (!check) {
    stop('Wrong format of minM and maxM. Please check!')
  }

  # Reshape minimal matrix
  res_element <- setdiff(rownames(maxM), rownames(minM))
  minM_row <-
    matrix(
      0,
      nrow = length(res_element),
      ncol = ncol(minM),
      dimnames = list(res_element, colnames(minM))
    )
  minM_col <-
    matrix(
      0,
      nrow = ncol(maxM),
      ncol = length(res_element),
      dimnames = list(c(rownames(minM), res_element), res_element)
    )
  minM <- cbind(rbind(minM, minM_row), minM_col)
  minM <- minM[rownames(maxM), rownames(maxM)]
  return(minM)
}
