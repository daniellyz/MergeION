#' Screening of LC-MS/MS chromatograms for metabolic features
#'
#' The function performs a fast non-targeted metabolic feature search in all input datafiles. It provides a metadata matrix for library_generator.
#'
#' @param raw_data_files A character vector of file names of chromatograms fOr feature screening. All files must have be in centroid-mode with mzML or mzMXL extension!
#' @param polarity Character 'positive' or 'negative'
#' @param rt_search Retention time search tolerance (in second) for detecting common metabolic features
#' @param ppm_search m/z search tolerance (in ppm) for detecting common metabolic features
#' @param baseline Numeric. Minimal absolute intensity to be considered as a metabolic feature.
#' @return
#' ref: Matrix object with 5 columns. Metadata containing all detected features
#'
#' @importFrom xcms xcmsSet group.nearest peakTable

MS1_screener <- function(raw_data_files, polarity = 'positive', rt_search = 12, ppm_search = 10, baseline = 1000){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  # Peak detection from MS1 scans of raw datafiles

  xset_full<- xcmsSet(raw_data_files, method="matchedFilter", mslevel = 1, polarity = polarity)
  mz_search = ppm_search/1000000*500
  xset_full <- group.nearest(xset_full, mzVsRTbalance=10, mzCheck=mz_search, rtCheck=rt_search)
  peak_info <- peakTable(xset_full,filebase="peakList")
  peak_info <- as.data.frame(peak_info)

  # Filter noise:
  NC = ncol(peak_info)
  I_info =peak_info[,(NC+1-length(raw_data_files)):NC]
  I_info[is.na(I_info)]=0
  I_info = apply(I_info,1,mean)
  valid = which(I_info>=baseline*5)
  peak_info = peak_info[valid,]

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
