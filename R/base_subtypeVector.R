
#' @rdname subtypeVector
#' @title subtypeVector
#' @description Automatic clustering based on expression matrix and gene expression profiles
#' @param geneSet List. Gene expression profiles - Co-expressed modules from WGCNA, etc
#' @param outlier.filter Integer. Set a small integer to remove outlier samples
#' @param k the number of clusters you want. Default is the same as the number of \code{geneSet} types
#' @inheritParams PAD
#' @importFrom luckyBase mergeMatrixDup
#' @importFrom luckyBase convert
#' @importFrom plyr ldply
#' @importFrom stats hclust
#' @importFrom stats dist
#' @importFrom stats cutree
#' @importFrom stringr str_extract
#' @importFrom WGCNA labels2colors
#' @import ComplexHeatmap
#' @return A list: \cr
#' \itemize{
#'   \item \code{Repeat} Basic parameters
#'   \item \code{Data} mean expression of each GEP in each sample cluster
#'   \item \code{logloss} negative log-likelihood function
#'   \item \code{MissGene} Missing genes
#'   \item \code{Plot} a \code{\link[ComplexHeatmap]{Heatmap}} object
#' }
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
subtypeVector <- function(expr,
                          geneSet,
                          cluster.method = c("ward.D2", 'complete')[1],
                          extra.annot = NULL,
                          # extra.annot = ComplexHeatmap::HeatmapAnnotation(),
                          plot.title = NULL,
                          k = NULL,
                          outlier.filter = 2,
                          verbose = T) {
  # Test
  if (F) {
    library(luckyBase)
    np <-
      c('tidyr', 'plyr', 'dplyr', 'stringr', 'ComplexHeatmap')
    Plus.library(np)

    # set.seed(200); expr <- readRDS("E:/Sync/@Analysis/PanCan_Data/Other/Level_02/tcga_RSEM_gene_tpm.rds"); expr <- expr[sample(1:ncol(expr),1000, replace = F)];
    expr <-
      readRDS("E:/Sync/@Analysis/PanCan_Data/Other/Level_02/tcga_RSEM_gene_tpm.rds")
    geneSet <-
      readRDS('E:/RCloud/database/Signature/report/Signature_PanCanWGCNA-Top10.rds')
    cluster.method = "ward.D2"
    k = NULL
    extra.annot = NULL
    plot.title = NULL
    verbose = T
    outlier.filter = 2
  }

  # Gene expression profiles (geneSet)
  if (T) {
    if (verbose)
      LuckyVerbose('subtypeVector: Data preparation...')
    geneSet_Gene <- unlist(geneSet)
    coGene <- intersect(geneSet_Gene, rownames(expr))
    # ID match test
    if (length(coGene) < length(geneSet_Gene)) {
      lack <- setdiff(geneSet_Gene, coGene)
      if (verbose)
        LuckyVerbose(
          'subtypeVector: Gene match:',
          sprintf('%.1f', length(coGene) * 100 / length(geneSet_Gene)),
          '%. Lack of ',
          paste0(lack, collapse = ', '),
          '.'
        )
    } else {
      lack <- NULL
      if (verbose)
        LuckyVerbose('subtypeVector: Gene match: 100%.')
    }
    for (i in 1:length(geneSet)) {
      # i=1
      geneSet[[i]] <- geneSet[[i]][geneSet[[i]] %in% coGene]
    }
    x <- as.matrix(expr)[as.character(unlist(geneSet)),]
    xZ <- t(scale(t(x), center = T, scale = T))
  }

  # Subtypes
  if (T) {
    if (verbose)
      LuckyVerbose('subtypeVector: Clustering...')
    dis <- dist(t(xZ), method = "euclidean")
    hc <- hclust(dis, method =  cluster.method)
    if (is.null(k)) {
      k <- length(geneSet)
    }
    subtype_vector <- cutree(hc, k)
  }

  # Alignment
  if (T) {
    if (verbose)
      LuckyVerbose('subtypeVector: Alignment...')
    gene_annot <- ldply(geneSet, cbind)
    colnames(gene_annot) <- c('module', 'gene')


    xZ_meanExpress <- mergeMatrixDup(
      x = xZ,
      mergeCol = T,
      fun_col = function(x)
        median(x, na.rm = T),
      refCol = as.character(subtype_vector[colnames(x)]),
      mergeRow = T,
      fun_row = function(x)
        median(x, na.rm = T),
      refRow = convert(rownames(x), 'gene', 'module', gene_annot),
      verbose = F
    )

    xZ_align <-
      apply(xZ_meanExpress, 2, function(x) {
        names(x)[x == max(x, na.rm = T)][1]
      })

    df_align <- data.frame(
      clusterSubtype = names(xZ_align),
      # geneSetubtype = str_extract(as.character(xZ_align), '[0-9]{1,4}'),
      geneSetubtype = as.character(xZ_align),
      stringsAsFactors = F
    )
    subtype_vector_aligned <-
      convert(subtype_vector,
              'clusterSubtype',
              'geneSetubtype',
              df_align)
    names(subtype_vector_aligned) <- names(subtype_vector)

    # Remove outliers
    subtype_vector_aligned <- subtype_vector_aligned[]
    target_cluster <-
      names(table(subtype_vector_aligned))[table(subtype_vector_aligned) > outlier.filter]
    abandon_cluster <-
      names(table(subtype_vector_aligned))[table(subtype_vector_aligned) <= outlier.filter]
    if (verbose)
      LuckyVerbose(
        'subtypeVector: Target clusters: ',
        paste(target_cluster, collapse = ', '),
        '. Abandoned cluster: ',
        paste(abandon_cluster, collapse = ', '),
        '.'
      )
    subtype_vector_aligned <-
      subtype_vector_aligned[subtype_vector_aligned %in% target_cluster]

  }

  # Heatmap
  if (T) {
    if (verbose)
      LuckyVerbose('subtypeVector: Heatmap...')

    # data
    # module_number <- str_extract(gene_annot$module, '[0-9]{1,4}')
    geneSet_vector <- as.character(gene_annot$module)
    geneSet_unique <- unique(geneSet_vector)
    geneSet_annot <- data.frame(
      geneSet = geneSet_unique,
      Color = c(mycolor, setdiff(scales::hue_pal()(
        length(geneSet_unique)
      ), mycolor))[1:length(geneSet_unique)],
      stringsAsFactors = F
    )
    xZ <-
      xZ[as.character(gene_annot$gene), names(subtype_vector_aligned)]


    # row annotation
    ra <- rowAnnotation(
      geneSet = geneSet_vector,
      col = list(geneSet = {
        col <-
          convert(geneSet_vector, 'geneSet', 'Color', geneSet_annot)
        names(col) <- geneSet_vector
        col
      }),
      annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
      annotation_name_side = "top"
    )

    ht_opt$message = FALSE # Turn off warning message
    suppressMessages(
      p <- Heatmap(
        xZ,
        name = 'z-score',
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_split = as.character(subtype_vector_aligned),
        row_split = geneSet_vector,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        column_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 8, fontface = "bold"),
        row_title = plot.title,
        row_title_side = 'left',
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        left_annotation = ra
      )
    )

    if (!is.null(extra.annot))
      p <- p %v% extra.annot
    # if(verbose) print(p)
  }

  # Output
  res <- list(
    Repeat = list(
      geneSet = geneSet,
      cluster.method = cluster.method,
      plot.title = plot.title,
      extra.annot = extra.annot
    ),
    Data = list(
      # meanExpressMatrix = xZ_meanExpress,
      expr = x[, names(subtype_vector_aligned)],
      subtype = translate_subtype(subtype_vector_aligned)
    ),
    MissGene = lack,
    Plot = p
  )
  if (verbose)
    LuckyVerbose('Done!')
  return(res)

}


####%%%%%%%%%%%%%%%%Assistant function%%%%%%%%%%%%%%%%####

#' @importFrom luckyBase convert
translate_subtype <- function(subtype_vector_aligned) {
  n <- names(subtype_vector_aligned)

  if (sum(grepl('Module_', subtype_vector_aligned)) == length(subtype_vector_aligned)) {
    # subtype = 'Module_xx'
    subtype_vector_aligned <-
      gsub('Module_', '', subtype_vector_aligned)
  } else {
    subtype_annot <- data.frame(
      raw = unique(subtype_vector_aligned),
      new = 1:length(unique(subtype_vector_aligned)),
      stringsAsFactors = F
    )
    subtype_vector_aligned <-
      convert(subtype_vector_aligned, 'raw', 'new', subtype_annot)
  }

  names(subtype_vector_aligned) <- n

  return(subtype_vector_aligned)

}
