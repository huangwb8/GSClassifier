

#' @description classify a vector by specified classifier
#' @param vector a vector
#' @param classifier  a list of classifier
#' @param cover if \code{cover=T}(default), informations not provided by \code{classifier} would be annotated as \code{Not Available}. Otherwise, they would be kept as raw.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## example
#' vector <- 1:100
#'
#' ## classifier exist
#' classifier1 <- list(
#'   lower = c(1:49),
#'   upper = c(50:90)
#' )
#' v1 <- classify(vector,classifier1)
#' table(v1)
#'
#' ## upper not exist
#' classifier2 <- list(
#'   lower = c(1:49),
#'   upper = c(101:900)
#' )
#' v1 <- classify(vector,classifier2)
#' table(v1)
classify <- function(vector,classifier,cover=T){
  vector <- as.character(vector)

  ## standard vector
  if(cover){
    vector1 <- rep("Not Available",length(vector))
  } else {
    vector1 <- vector
  }

  ## replace
  for(i in 1:length(classifier)){ # i = 1
    n.i <- names(classifier)[i]
    c.i <- as.character(classifier[[i]])
    vector1[vector %in% c.i] <- n.i
  }
  return(vector1)
}



#' @title PAD
#' @description Unsupervised pan-immune activation/dysfunction (PAD) subtypes of
#'   gastric cancer (or other solid tumor) sample based on RNA-Seq/microarray
#'   data
#' @param expr RNA expression matrix. Samples in col and ENSEMBL genes in row.
#' @param PIAM ID of pan-immune activation genes. IF \code{NULL}, use default
#'   gene set.
#' @param PIDG ID of pan-immune dysfunction genes. IF \code{NULL}, use default
#'   gene set.
#' @param cluster.method One of \code{'ward.D2'}, \code{'complete'} or other
#'   methods in \code{\link[stats]{hclust}}. If \code{'randomForest'} was set,
#'   the \code{\link[randomForest]{randomForest}} function would used to
#'   classify the samples.
#' @param rF.para The parameters in \code{\link[randomForest]{randomForest}}
#'   function.
#' @param plot.title The title of heatmap report.
#' @param verbose Whether to show heatmap in the process.
#' @param extra.annot Extra top annotation. The same order as colnames of \code{expr}.
#' @inheritParams callEnsemble
#' @importFrom stats dist hclust cutree
#' @importFrom randomForest randomForest
#' @importFrom cluster pam
#' @importFrom dplyr left_join arrange
#' @importFrom tidyr %>%
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation Heatmap
#' @importFrom grid gpar
#' @details This function is used for unsupervised classification of raw data,
#'   which is pivotal for the following supervised machine learning. Empirically, the \code{'ward.D2'} method could be useful and
#'   high-speed for simple gene signatrues (like PAD classifier). Random forest
#'   is a powerful stragety and may act well in larger dataset or complex gene
#'   signatures.
#' @examples
#' extra.annot = HeatmapAnnotation(
#' Dataset = id.dataset,
#' col = list(
#'   Dataset = if(T){
#'     l <- mycolor[1:length(unique(id.dataset))];
#'     names(l) <- unique(id.dataset);
#'     l}
#' ),
#' annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
#' show_legend = T
#' )
#'
#' res1 <- PAD(
#'   expr = dm.combat.tumor,
#'   PIAM = piam,
#'   PIDG = pidg,
#'   plot.title = 'PanSTAD',
#'   cluster.method = 'ward.D2',
#'   subtype = 'PAD.train_20200110',
#'   extra.annot = extra.annot,
#'   verbose = T
#' )
#'
#' # randomForest: time-consuming in large cohorts
#' res2 <- PAD(
#'   expr = dm.combat.tumor,
#'   PIAM = piam,
#'   PIDG = pidg,
#'   cluster.method = 'randomForest',
#'   rF.para = list(
#'     seed = c(2020,485,58,152),
#'     ntree = c(1000,1000),
#'     k=c(2,2)
#'   ),
#'   subtype = 'PAD.train_20200110',
#'   extra.annot = extra.annot,
#'   plot.title = 'PanSTAD',
#'   verbose = T
#' )
#' @export
PAD <- function(
  expr,
  PIAM = NULL,
  PIDG = NULL,
  cluster.method = c("ward.D2",
                     'complete',
                     'randomForest')[1],
  rF.para = list(
    seed = c(2020,485,4,8),
    ntree = c(300,300),
    k=c(2,2)
  ),
  extra.annot = NULL,
  # extra.annot = ComplexHeatmap::HeatmapAnnotation(),
  plot.title = NULL,
  subtype = 'PAD.train_20200110',
  verbose = T
){

  ## Test
  if(F){
    expr <- t(dataExpr)
    PIAM = piam
    PIDG = pidg
    cluster.method = c("ward.D2",'complete','randomForest')[1]
    rF.para = list(
      seed = c(2020,485,4,8),
      ntree = c(300,300),
      k=c(2,2)
    )
    plot.title = 'TCGA-STAD'
    verbose = T
  }

  ## Package
  # library(lucky)
  # nd <- c('NbClust','dplyr','ComplexHeatmap','reshape2','randomForest','cluster'); lucky::Plus.library(nd)

  ## PIAM & PIDG
  l <- readRDS(system.file("extdata", paste0(subtype, '.rds'), package = "GSClassifier"))
  if(is.null(PIAM)) {
    PIAM <- l$geneSet$PIAM
    cat('Use default PIAM...','\n')
  }
  if(is.null(PIDG)){
    PIDG <- l$geneSet$PIDG
    cat('Use default PIDG...','\n')
  }

  ## Data
  coGene <- intersect(c(PIAM,PIDG),rownames(expr))
  # ID match test
  if(length(coGene) < length(c(PIAM,PIDG))){
    lack <- setdiff(c(PIAM,PIDG),coGene)
    cat('Gene match:',sprintf('%.1f',length(coGene)*100/length(c(PIAM,PIDG))),'%. Lack of ',paste0(lack,collapse = ', '),'.','\n')
  } else {
    lack <- NULL
    cat('Gene match: 100%.','\n')
  }
  PIDG <- PIDG[PIDG %in% coGene]
  PIAM <- PIAM[PIAM %in% coGene]
  x <- as.matrix(expr)[c(PIAM,PIDG),]
  xZ <- t(scale(t(x),center = T,scale = T))

  ## PIAM subtypes
  if(T){
    xZ.piam <- as.data.frame(t(xZ[PIAM,]))
    if(cluster.method != 'randomForest'){
      # Use non-randomForest strategy
      # numComplete <- NbClust(xZ.piam,
      #                        distance = "euclidean",
      #                        min.nc = 2,
      #                        max.nc= 4,
      #                        method = "complete",
      #                        index = "all")
      dis <- dist(xZ.piam, method = "euclidean")
      hc <- hclust(dis, method =  cluster.method) # "ward.D2"
      comp2 <- cutree(hc, 2)
    } else {
      # Use randomForest strategy
      set.seed(rF.para$seed[1]); rf <- randomForest(x = xZ.piam, ntree = rF.para$ntree[1], proximity = T)
      # dim(rf$proximity)
      # View(rf$proximity[1:5, 1:5])
      # importance(rf)
      dissMat <- sqrt(1 - rf$proximity)
      set.seed(rF.para$seed[2]); pamRF <- pam(dissMat, k = rF.para$k[1])
      comp2 <- pamRF$clustering
    }
    p1 <- sum(colMeans(xZ.piam[names(comp2)[comp2 == 1],]))
    p2 <- sum(colMeans(xZ.piam[names(comp2)[comp2 == 2],]))

    if(p1>p2){
      id_high <- names(comp2)[comp2 == 1]
      h <- rep('piam_high',length(id_high));names(h) <- id_high
      id_low <- names(comp2)[comp2 == 2]
      l <- rep('piam_low',length(id_low));names(l) <- id_low
    } else {
      id_high <- names(comp2)[comp2 == 2]
      h <- rep('piam_high',length(id_high));names(h) <- id_high
      id_low <- names(comp2)[comp2 == 1]
      l <- rep('piam_low',length(id_low));names(l) <- id_low
    }

    hl <- c(h,l)

    testPIAM <- data.frame(
      ID = names(hl),
      PIAM = hl,
      stringsAsFactors = F
    )
  }

  ## PIDG subtypes
  if(T){
    xZ.PIDG <- as.data.frame(t(xZ[PIDG,]))
    if(cluster.method != 'randomForest'){
      # Use non-randomForest strategy
      dis <- dist(xZ.PIDG, method = "euclidean")
      hc <- hclust(dis, method =  cluster.method) # "ward.D2"
      comp2 <- cutree(hc, 2)
    } else {
      # Use randomForest strategy
      set.seed(rF.para$seed[3]); rf <- randomForest(x = xZ.PIDG, ntree = rF.para$ntree[2], proximity = T)
      dissMat <- sqrt(1 - rf$proximity)
      set.seed(rF.para$seed[4]); pamRF <- pam(dissMat, k = rF.para$k[2])
      comp2 <- pamRF$clustering
    }

    # Get high/low subgroup
    p1 <- sum(colMeans(xZ.PIDG[names(comp2)[comp2 == 1],]))
    p2 <- sum(colMeans(xZ.PIDG[names(comp2)[comp2 == 2],]))

    if(p1>p2){
      id_high <- names(comp2)[comp2 == 1]
      h <- rep('pidg_high',length(id_high));names(h) <- id_high
      id_low <- names(comp2)[comp2 == 2]
      l <- rep('pidg_low',length(id_low));names(l) <- id_low
    } else {
      id_high <- names(comp2)[comp2 == 2]
      h <- rep('pidg_high',length(id_high));names(h) <- id_high
      id_low <- names(comp2)[comp2 == 1]
      l <- rep('pidg_low',length(id_low));names(l) <- id_low
    }

    hl <- c(h,l)

    testPIDG <- data.frame(
      ID = names(hl),
      PIDG = hl,
      stringsAsFactors = F
    )
  }

  ## PIAM+PIDG
  df <- left_join(testPIAM,testPIDG,by='ID')
  a <- paste(df$PIAM,df$PIDG,sep = '-') # table(df$`PAD subtype`)
  df$`PAD subtype` <- ifelse(
    a == 'piam_high-pidg_high','PAD-I',ifelse(
      a == 'piam_high-pidg_low','PAD-II',ifelse(
        a == 'piam_low-pidg_high','PAD-III',
        'PAD-IV')))
  df$PIAM <- gsub('piam_','',df$PIAM)
  df$PIDG <- gsub('pidg_','',df$PIDG)
  # df <- arrange(df,`PAD subtype`) %>% as.data.frame()
  df <- df[match(colnames(xZ),df$ID),] # aligned to expr

  ## TMEscore-like: Activated immune score
  # This score is not well for estimation of gene expression.
  if(F){
    # PC1 of  PIAM
    xZ_PIAM <- xZ[PIAM,]
    pca.mt <- prcomp(xZ_PIAM,center = F,scale. = F)
    PC1_PIAM <- pca.mt[["rotation"]][,1]

    # PC1 of PIDG
    xZ_PIDG <- xZ[PIDG,]
    pca.mt <- prcomp(xZ_PIDG,center = F,scale. = F)
    PC1_PIDG <- pca.mt[["rotation"]][,1]

    # Match
    PC1_PIAM <- PC1_PIAM[names(PC1_PIDG)]
    df_score <- data.frame(
      ID = names(PC1_PIAM),
      AIscore = PC1_PIAM - PC1_PIDG, # Activated immune score
      stringsAsFactors = F
    )

    df2 <- left_join(df,df_score,by='ID') %>%
      arrange(`PAD subtype`,desc(AIscore)) %>%
      as.data.frame()
  }

  ## Heatmap
  if(T){
    rownames(df) <- as.character(df$ID)
    # xZ <- xZ[,rownames(df)]
    # df <- df[match(table = as.character(df$ID),colnames(xZ)),]

    # col annotation
    ha = HeatmapAnnotation(
      PIAM = df$PIAM,
      PIDG = df$PIDG,
      `PADm subtype` = df$`PAD subtype`,
      col = list(
        PIAM = c('high' = 'red','low' = 'blue'),
        PIDG = c('high' = 'red','low' = 'blue'),
        `PADm subtype` = c(
          'PAD-I' = "#F8766D",
          'PAD-II' = "#7CAE00",
          'PAD-III' = "#00BFC4",
          'PAD-IV' = "#C77CFF"
        )
      ),
      annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
      show_legend = T
    )

    # row annotation
    ra <- rowAnnotation(
      Gene = classify(
        rownames(xZ),
        list(
          PIAM = PIAM,
          PIDG = PIDG
        ),
        cover = F
      ),
      col = list(
        Gene = c(
          PIAM = 'red',
          PIDG = 'blue'
        )
      ),
      annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
      annotation_name_side = "top"
    )

    p <- Heatmap(xZ,
                 name = 'z-score',
                 cluster_rows = T,
                 cluster_columns = F,
                 show_column_names = F,
                 show_row_names = F,
                 column_split = df$`PAD subtype`,
                 row_km = 2,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_columns = "euclidean",
                 column_names_gp = gpar(fontsize = 12, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                 row_title = plot.title,
                 row_title_side = 'left',
                 row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                 top_annotation = ha,
                 left_annotation = ra
    )
    if(!is.null(extra.annot)) p <- p %v% extra.annot
    if(verbose) print(p)

  }

  ## Output
  res <- list(
    Repeat = list(
      expr = x,
      PIAM = PIAM,
      PIDG = PIDG,
      cluster.method = cluster.method,
      rF.para = rF.para,
      plot.title = plot.title,
      subtype = subtype,
      extra.annot = extra.annot
    ),
    Data = df,
    MissGene = lack,
    Plot = p
  )
  cat('Done!','\n')
  return(res)
}









