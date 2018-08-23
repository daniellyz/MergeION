#' Screening of LC-MS/MS chromatograms for metabolic features
#'
#' The function performs a fast non-targeted metabolic feature search in all input datafiles. It provides a metadata matrix for library_generator.
#'
#' @param raw_data_files A character vector of file names of chromatograms fOr feature screening. All files must have be in centroid-mode with mzML or mzMXL extension!
#' @param polarity Character 'positive' or 'negative'
#'
#' @return
#' ref: Matrix object with 5 columns. Metadata containing all detected features
#'
#' @importFrom xcms xcmsSet group peaks

MS1_screener <- function(raw_data_files, polarity = 'positive'){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  # Peak detection from MS1 scans of raw datafiles

  xset_full<- xcmsSet(raw_data_files, method="matchedFilter", mslevel = 1, polarity = polarity)
  xset_full <- group(xset_full)
  peak_info <- peaks(xset_full)
  peak_info <- as.data.frame(peak_info)

  # Create metadata

  PI = nrow(peak_info)
  ref = data.frame(matrix(nrow =PI, ncol = 5))
  colnames(ref) = c("PEPMASS","RT","IONMODE","ADDUCT","ID")
  ref$PEPMASS = peak_info$mz
  ref$RT = peak_info$rt/60
  if (polarity == "positive"){
    ref$IONMODE = "Positive"
    ref$ADDUCT = "M+H"}
  if (polarity == "negative"){
    ref$IONMODE = "Negative"
    ref$ADDUCT = "M-H"}
  ref$ID = paste0("XCMS_", 1:PI)
  return(ref)
}
