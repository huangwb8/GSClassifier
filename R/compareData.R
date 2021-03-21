


#' @title Compare data from the same source but different scale
#' @description Compare data from the same source but different scale
#' @param X1 Matrix with sample col and gene row
#' @param X2 Matrix with sample col and gene row
#' @param names Label of X1 and X2
#' @param width,height Plot parameters
#' @param savePath Plot saved path
#' @importFrom plyr alply
#' @importFrom ggplot2 ggplot geom_point aes stat_smooth labs theme_bw theme
#'   element_text element_rect
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## PAD genes
#' PAD <- readRDS(system.file("extdata", 'PAD.rds', package = "GSClassifier"))
#' tgene <- as.character(unlist(PAD$geneSet))
#' X1 <- stad_fpkm[tgene,]
#' X2 <- stad_tpm[tgene,]
#'
#' ## Compare FPKM with TPM matrix
#' res <- compareData(X1,X2,
#'                    names=c('FPKM','TPM'),
#'                    width = 10,
#'                    height = 10,
#'                    savePath = '.')
#' res2 <- compareData(log2(X1+1),log2(X2+1),
#'                     names=c('log2FPKM','log2TPM'),
#'                     width = 10,
#'                     height = 10,
#'                     savePath = '.')
#' @export
compareData <- function(X1,
                        X2,
                        names=c('X1','X2'),
                        width = 10,
                        height = 10,
                        savePath = '.'){

  ## Data
  coGene <- intersect(rownames(X1),rownames(X2))
  coSample <- intersect(colnames(X1),colnames(X2))
  X1 <- X1[coGene,coSample]
  X2 <- X2[coGene,coSample]

  ## Correlation of expression for every gene
  compareOneGene <- function(X1,X2,gene){ # gene='ENSG00000136167'
    x1 <- as.numeric(X1[gene,])
    x2 <- as.numeric(X2[gene,])
    res <- cor.test(x1, x2,
                    alternative = "two.sided",
                    method = "spearman",
                    conf.level = 0.95)
    # res2 <- data.frame(
    #   Cor = res$estimate,
    #   P.value = res$p.value
    # )
    size = 30
    p <- ggplot(data.frame(x1=x1,x2=x2),aes(x=x1,y=x2)) +
      geom_point(alpha = 0.5,size = 11,colour = "#FB8072") +
      stat_smooth(formula = y ~ x,method = 'glm') +
      labs(x = names[1],
           y = names[2],
           title = gene) +
      theme_bw() +
      theme(
        plot.title = element_text(size = size,colour = "black",face = "bold",hjust = 0.5),
        axis.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        axis.title.x = element_text(size = size,colour = "black",face = "bold"),
        axis.title.y = element_text(size = size,colour = "black",face = "bold"),
        legend.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.title = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.position='right',
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = size/15*12,colour = "black",face = "bold")
      )

    # Output
    l <- list(
      Data = res,
      Plot = p
    )
    return(l)

  }
  l <- alply(as.matrix(coGene),1,function(x) compareOneGene(X1,X2,x))
  names(l) <- coGene

  ## Merge plot
  cairo_pdf(paste0(savePath,'/Plot of compareData_',paste0(names,collapse = '-'),'.pdf'),width = width,height = height,onefile = T)
  for(i in 1:length(l)) print(l[[i]]$Plot)
  dev.off()

  ## Output
  return(l)
}


















