convert.orthologs <- function(genes, 
                              table = NULL,  
                              from="Gene.HS", 
                              to="Gene.MM") {
  if (is.null(table)) {
    load("Hs2Mm.convert.table.RData")
  }
  genes.select <- genes[genes %in% table[,from]]
  
  #Convert
  ortho.genes <- table[,to][match(genes, table[,from])]
  return(ortho.genes)
}