writeMGF2 <- function(splist, metadata, con) {
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  con <- file(description = con, open = "at")
  on.exit(close(con))
  
  N=nrow(metadata)
  C=ncol(metadata)
  labels=colnames(metadata)
  for (i in 1:N) {
    .cat("\nBEGIN IONS\n")
    for (j in 1:C){
      .cat(labels[j],"=",metadata[i,j],"\n")}
    sp=splist[[i]]
    .cat(paste(sp[,1],"\t",sp[,2], collapse = "\n"))
    .cat("\nEND IONS\n")
  }
}
