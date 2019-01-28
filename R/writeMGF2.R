#' Writing metadata and a list of spectra to mgf files
#'
#' The function writes a dataframe of metadata and a list of spectra to a mgf file
#'
#' @param library List that contains 2 elements.
#' \itemize{
#'  \item{"sp" ~ List of spectra ready to be written. Each spectrum is a data matrix with two columns: m/z and intensity.}
#'  \item{"metadata" ~ Data frame containing metadata corresponding to each "scan" in splist.}
#'  }
#' @param con Name of the output library, the file extension must be mgf
#'
#' @examples
#' # Read a library file:
#' library = readMGF2("library_V2.mgf")
#' # Add new metadata "RESOLUTION = HIGH" to all scans:
#' library$metadata$RESOLUTION = "HIGH"
#' # Write into a new mgf file:
#' writeMGF2(library,"library_V2_bis.mgf")
#'
#' @importFrom tools file_ext
#' @export
#'
writeMGF2 <- function(library, con) {

 if (is.list(library)){
    if (length(library)!=2 || (!is.list(library$sp)) || !is.data.frame(library$metadata)){
      stop("Please make sure your input library is a valid output of library_generator() or readMGF2()!")
    }} else {stop("Please make sure your input library is a valid output of library_generator() or readMGF2()!")
  }

 if (is.character(con)){
    if (file_ext(con)!="mgf"){
      stop("The file extension of your input library must be mgf!")
    }} else {
      stop("The input must be the name of the mgf file!")}

  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }

  con <- file(description = con, open = "wt")
  on.exit(close(con))

  metadata = library$metadata
  splist = library$sp
  N=nrow(metadata)
  C=ncol(metadata)
  labels=colnames(metadata)
  for (i in 1:N) {
    .cat("\nBEGIN IONS\n")
    for (j in 1:C){
      .cat(labels[j],"=",as.character(metadata[i,j]),"\n")}
    sp=splist[[i]]
    .cat(paste(sp[,1],"\t",sp[,2], collapse = "\n"))
    .cat("\nEND IONS\n")
  }
}