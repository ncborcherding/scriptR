

prepare.data <- function(dir, 
                         organism = "human", 
                         normalization.method = "LogNormalize", 
                         use.sct = FALSE,
                         variable.features = 2000,
                         regress.variables = c("mito", "ribosomes", 'cell.cycle')) {
  
  seurat.list <- lapply(dir, functions(x) {
    Read10X(x)
  })
  
  if (organism == human) {
    mito.pattern <- "^MT-"
    ribo.pattern <- "^RBS|RPL"
    s.genes <- Seurat::cc.genes.updated.2019$s.genes
    g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  } else if (organism == "mouse") {
    mito.pattern <- tolower("^MT-")
    ribo.pattern <- paste0(str_to_title("^RBS") ,"|", str_to_title("RPL"))
    s.genes <- str_to_title(Seurat::cc.genes.updated.2019$s.genes)
    g2m.genes <- str_to_title(Seurat::cc.genes.updated.2019$g2m.genes)
  } else {
    stop("Currently only supporting human or mouse data preparation")
  }
  
  for(x in seq_along(seurat.list)) {
    seurat.list[[x]] <- CellCycleScoring(seurat.list[[x]], 
                                         s.features = s.genes, 
                                         g2m.features = g2m.genes, 
                                         set.ident = FALSE)
    seurat.list[[x]][["percent.mt"]] <- PercentageFeatureSet(seurat.list[[x]], pattern = mito.pattern)
    seurat.list[[x]][["percent.rb"]] <- PercentageFeatureSet(seurat.list[[x]], pattern = ribo.pattern)
  }

  if (!use.sct) {
    seurat.list <- lapply(X = seurat.list, FUN = function(x) {
      x <- NormalizeData(x, 
                         normalization.method = normalization.method)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = variable.features)
    })
  } else {
    seurat.list <- lapply(X = ifnb.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = variable.features)
    seurat.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
  }

}