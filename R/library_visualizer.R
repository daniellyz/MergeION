#' Visualize selected mass spectra in the spectral library
#'
#' The function plots all mass spectra found by library query
#'
#' @param query A list of four objects SELECTED_SCANS, LEFT_SCANS, SELECTED_LIBRARY and LEFT_LIBRARY. It should be the output of library_manager().
#' @param query.mode A character. Either "found" or "remained".
#' @param max.plot An integer. At least 1, maximal number of plots per page
#'
#' @examples
#' # Search library using query command lines:
#' query = library_manager(library1,query=c("IONMODE=Positive","MSLEVEL=2","RT=1.2"), logical="AND", rt_search=6)
#'
#' # Visualize scans found:
#' library_visualizer(query, query.mode="found")
#'
#' @export
#'
#' @importFrom MSnbase fData readMgfData
#' @importFrom graphics plot title legend abline text
#'
library_visualizer<-function(query, query.mode=c("found","remained"), max.plot=0){

  options(stringsAsFactors = FALSE)
  options(warn=-1)
  max_display = 5 # Display text of 5 most abundant mass peaks and higher than 5%

  #################
  ### Check inputs:
  #################

  if (missing(query)){
    stop("Please provide the output of library_manager() as input!")}

  if (is.list(query)){
    if (length(query)!=4 || (!is.list(query$SELECTED_LIBRARY)) || (!is.list(query$LEFT_LIBRARY))){
      stop("Please make sure your input library is a valid output of library_manager()!")
    }}

  query.mode = match.arg(query.mode,choices=c("found","remained"),several.ok = FALSE)

  ##################################
  ### Reading from spectral library:
  ##################################

  if (query.mode == "found"){
    metadata = query$SELECTED_LIBRARY$metadata
    spectrum_list = query$SELECTED_LIBRARY$sp
    if (max.plot==0){max.plot = nrow(metadata)}
  }

  if (query.mode == "remained"){
    metadata = query$LEFT_LIBRARY$metadata
    spectrum_list = query$LEFT_LIBRARY$sp
    if (max.plot==0){max.plot = nrow(metadata)}
  }

  ##############
  ### Visualize:
  ##############

  NI = nrow(metadata)
  par(mfrow=c(max.plot,1))
  if (NI>0){
    for (i in 1:NI){
      spectrum = spectrum_list[[i]]
      prec_mz = round(as.numeric(metadata$PEPMASS[i]),3)

      if (metadata$MSLEVEL[i]==1){
        xrange = c(prec_mz-2,prec_mz+8)
        ranges = which((spectrum[,1]>=prec_mz-2) & (spectrum[,1]<=prec_mz+8))
        yrange = c(0, max(spectrum[ranges,2])*1.2)}
      if (metadata$MSLEVEL[i]==2){
        xrange = c(50,prec_mz+8)
        ranges = which((spectrum[,1]>=50) & (spectrum[,1]<=prec_mz+8))
        yrange = c(0, max(spectrum[ranges,2])*1.2)}

      spectrum = spectrum[ranges,]
      if (length(ranges)==1){spectrum = matrix(spectrum, ncol=2)}

      plot(spectrum[,1],spectrum[,2], type = "h", xlim = xrange, ylim = yrange,
           xlab = "m/z", ylab = "Intensity", font.lab=2)

      kkk = min(nrow(spectrum), max_display)
      max_pics = order(spectrum[,2], decreasing = T)[1:kkk]
      max_pics = intersect(max_pics,which(spectrum[,2]>=max(spectrum[,2])*0.05))
      text_pics = spectrum[max_pics,1]
      int_pics = spectrum[max_pics,2]*1.1
      text(text_pics, int_pics, as.character(round(text_pics,3)))

      title(paste0("ID: ",metadata$ID[i]))
      legend("topright", bty = "n",
            legend = paste0(
            "Precursor: ", prec_mz,  "\n",
           #    "RT (min): ", round(as.numeric(metadata$RT[ind]),2),  "\n",
           #    "ADDUCT: ", metadata$ADDUCT[ind], "\n",
            "MS Level: ", metadata$MSLEVEL[i]))

     if (metadata$MSLEVEL[i]==2){abline(v = prec_mz, col = "red")}
    }} else {
     print("No scans available for ID specified!")
    }
  }
