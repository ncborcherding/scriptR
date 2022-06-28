#' Calculate gene set enrichment scores between variables 
#'
#' Function is meant to cluster/sample/variable of interest
#' most closely related to in terms of gene enrichment
#' @param marker.list List of Differential Genes
#' @param sc Single-Cell Object
#' @param baseline Variable to compare
#' @param top.genes Number of genes used to make rank lists
#' @param return.table Return the results as a table
#' 
#' @importFrom dplyr bind_rows
#' @import ggplot2
#'
#' 
#' @examples 
#' DEG.list <- diffgene.lister(Macs2, 
#'                             ident.use = "idents", 
#'                             comparing.variable = "sample",
#'                             compare.by = c("Bleo", "PBS"))
#'                             
#' differential.enrich.results <- (DEG.list, 
#'                                  Macs2,
#'                                  baseline = "IntMacs1", 
#'                                  top.genes = 200)
#' 
#' @export
#' @author Nick Borcherding
#'
#' @return graph or table of the enrichment results
differential.enrich.results <- function(marker.list, 
                          sc = NULL, 
                          baseline = NULL, 
                          top.genes = 200, 
                          return.table = FALSE) {
  which.baseline <- which(names(marker.list) == baseline)
  gene.markers <- marker.list[[which.baseline]]
  
  gene.markers <- gene.markers[gene.markers$p_val_adj < 0.05,]
  gene.markers$trend <- ifelse(gene.markers[,"avg_log2FC"] > 0, "Up", "Down")
  gene.markers <- gene.markers[order(gene.markers$avg_log2FC, decreasing = TRUE),]
  gene.list <- split(gene.markers, gene.markers[,"trend"])
  for (i in seq_along(gene.list)) {
    gene.list[[i]] <- rownames(gene.list[[i]])
  }
  
  compare.list <- marker.list
  compare.list[which.baseline] <- NULL
  #####Replace processEnrich here
  out.list <- list()
  for(j in seq_along(compare.list)) {
    comp <- compare.list[[j]]
    up <- processEnrich(comparison = comp, sc = sc, genes = gene.list[[1]], baseline = baseline, top.genes)
    down <- processEnrich(comparison = comp, sc = sc, genes = gene.list[[2]], baseline = baseline, top.genes)
    out <- rbind.data.frame(up,down)
    out$Direction <- c("up", "down")
    out.list[[j]] <- out
    
  }
  out <- bind_rows(out.list)
  out$Description <- rep(names(compare.list), each = 2)
  out$Description <- paste0(out$Description, "_", out$Direction)
  out <- out[which(!is.na(out$ID)),]
  out$p.adjust <- p.adjust(out$pvalue)
  out$fill <- ifelse(out$p.adjust<0.05,"p.adj < 0.05",
                     ifelse(out$p.adjust>=0.05 & out$qvalue < 0.05,"q.value < 0.05","pvalue < 0.01"))
  if(return.table) {
    return(out)
  }
  plot <- ggplot(out, aes(reorder(Description, GeneRatio), GeneRatio, fill = fill)) +
    geom_col(width = .5) +
    scale_fill_manual(values=c("pvalue < 0.01" = "#F8766D", 
                               "p.adj < 0.05" = "#619CFF", "q.value < 0.05" = "#00BA38")) +
    coord_flip() +
    labs(x="Pathway", y="GeneRatio",
         title=baseline) + 
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "grey", size=1) +
    theme_minimal()
  return(plot)
}



#' Internal function to orgnaize enrichment results
#' @importFrom dplyr slice mutate arrange %>%
#' @importFrom clusterProfiler enricher
processEnrich <- function(comparison, 
                          sc,
                          genes, 
                          baseline, 
                          top.genes = top.genes) {
  require(clusterProfiler)
  universe <- rownames(sc@assays$RNA@counts)
  if (length(genes) > top.genes) {
    gen <- genes[seq_len(top.genes)]
  } else {
    gene <- genes
  }
  t2g <- data.frame(cellName = baseline, geneID = gen)
  x <- data.frame(cellName = "no", geneID = universe)
  t2g <- rbind(t2g,x)
  t2g <- na.omit(t2g)
  
  compare.tmp <- comparison[comparison$p_val_adj < 0.05,]
  compare.tmp <- compare.tmp[compare.tmp$avg_log2FC > 0,]
  compare.genes <- data.frame(geneID = rownames(compare.tmp), FC = compare.tmp$avg_log2FC)
  geneList1 = compare.genes[,2]
  names(geneList1) = as.character(compare.genes[,1])
  geneList1 = sort(geneList1, decreasing = TRUE)
  deg <- names(geneList1)
  if (length(intersect(deg, gen)) == 0) {
    tmp <- data.frame(matrix(ncol = 12, NA))
    colnames(tmp) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "fill","GeneRatio1", "GeneRatio2") 
    return(tmp)
  } else {                                    
    tmp <- enricher(gene     = deg,
                    TERM2GENE = t2g,
                    universe = universe,
                    maxGSSize = 1000)
    tmp <-  as.data.frame(tmp@result)
    rownames(tmp) <- c()
    tmp <-  tmp %>% 
      arrange(qvalue) %>%
      slice(c(1:15)) %>%
      mutate(fill = ifelse(p.adjust<0.05,"p.adj < 0.05",
                           ifelse(p.adjust>=0.05 & qvalue < 0.05,"q.value < 0.05","pvalue < 0.01")))
    
    tmp$GeneRatio1 <- stringr::str_extract(tmp$GeneRatio, "(.*)/")
    tmp$GeneRatio1 <- gsub("\\/", "", tmp$GeneRatio1)
    tmp$GeneRatio1 <- as.numeric(tmp$GeneRatio1)
    tmp$GeneRatio2 <- stringr::str_extract(tmp$GeneRatio, "/(.*)")
    tmp$GeneRatio2 <- gsub("\\/", "", tmp$GeneRatio2)
    tmp$GeneRatio2 <- as.numeric(tmp$GeneRatio2)
    tmp$GeneRatio <- tmp$GeneRatio1 / tmp$GeneRatio2
    tmp$GeneRatio <- tmp$GeneRatio * 100
  }
  return(tmp)
}

#' Generate a list of differential genes
#'
#' List of MAST-based differential genes to use for further analysis.
#' 
#' @param sc Single-Cell Object
#' @param idents.use Variable to set single-cell identity to use
#' @param comparing.variable meta variable to use for comparison
#' @param compare.by levels of the comparing variable to use
#' @param latent.vars latent variables to use for the MAST model

#' @importFrom dplyr bind_rows
#' @import ggplot2
#'
#' 
#' @examples 
#' DEG.list <- diffgene.lister(Macs2, 
#'                             ident.use = "idents", 
#'                             comparing.variable = "sample",
#'                             compare.by = c("Bleo", "PBS"))
#'                             
#' @export
#' @author Nick Borcherding
#'
#' @return graph or table of the enrichment results
diffgene.lister <- function(sc, 
                            ident.use = "idents", 
                            comparing.variable = "sample",
                            compare.by = c("Bleo", "PBS"), 
                            latent.vars = NULL) {
  require(Seurat)
  marker.list <- list()
  if (ident.use != "idents") {
    Idents(sc) <- ident.use
  }
  cells <- rownames(sc[[]])[which(sc@meta.data[,comparing.variable] %in% compare.by)]
  subset.tmp <- subset(sc, cells = cells)
  subset.tmp <- SplitObject(sc, split.by = "ident")
  
  for (i in seq_along(subset.tmp)) {
    if (is.null(latent.vars)) {
      mark <- FindMarkers(subset.tmp[[i]], 
                          group.by = comparing.variable, 
                          ident.1 = compare.by[1], 
                          ident.2 = compare.by[2], 
                          test.use = "MAST", 
                          pseudocount.use = 0.1)
    } else {
      mark <- FindMarkers(subset.tmp[[i]], 
                          group.by = comparing.variable, 
                          ident.1 = compare.by[1], 
                          ident.2 = compare.by[2], 
                          test.use = "MAST", 
                          pseudocount.use = 0.1, 
                          latent.vars = latent.vars)
    }
    marker.list[[i]] <- mark
    rm(mark)
  }
  names(marker.list) <- names(subset.tmp)
  return(marker.list)
}