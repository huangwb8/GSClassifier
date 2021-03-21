

#' @title modelData
#' @description Training and internal validation data selection
#' @param design A design object with ID rownames
#' @param id.col Colname of ID
#' @param variable Split factor for balanced selection
#' @param Prop the proportion of training cohort
#' @param seed random seed
#' @importFrom plyr ddply dlply
#' @export
modelData <- function(design,
                      id.col = 'ID',
                      variable = c('platform','PAD subtype'),
                      Prop = 0.6,
                      seed = 2020){

  ## Test
  if(F){
    design <- design3 # x <- expr
    variable = c('Dataset','PAD_subtype')
    Prop = 0.6
    seed = 2020
  }

  ## Data
  a <- ddply(design,
             .variables = variable,
             plyr::summarize,
             Size = length(platform))
  set.seed(2020); seeds <- sample(1:10000,nrow(a),replace = F)


  ## Subtset function
  getOneData <- function(x,Prop = 0.6,seed.i=seeds[1]){

    ## Test
    if(F){
      x <- expr[,design$ID[design$Dataset == 'GSE84437' & design$PAD_subtype == 'PAD-IV']]
    }

    ## Sample a subset
    set.seed(seed.i); idx <- sample(1:nrow(x), floor(nrow(x) * Prop) + 1, replace = F)

    ## Matrix
    if(length(idx) == 1){ # idx=1
      x2 <- x[idx,]
      x3 <- matrix(x2,nrow = 1,byrow = F,dimnames = list(rownames(x)[idx],names(x2)))
    } else {
      x3 <- x[idx,]
    }

    ## Output
    return(x3)

  }

  ## Training dataset
  l <- dlply(design,.variables = variable)
  L <- list()
  for(i in 1:length(l)){ # i=1
    x.i <- l[[i]]
    L[[names(l)[i]]] <- getOneData(x.i,Prop = Prop,seed.i=seeds[i])
  }
  train <- do.call('rbind',L)
  rownames(train) <- as.character(train[,id.col])

  ## Validation dataset
  valid <- design[setdiff(rownames(design),rownames(train)),]

  ## Output
  res <- list(
    Repeat = list(
      design = design,
      id.col = id.col,
      variable = variable,
      Prop = Prop,
      seed =  seed
    ),
    Data = list(
      Train = train,
      Valid = valid
    )
  )
  return(res)
}
