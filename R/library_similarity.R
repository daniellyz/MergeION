#' Matching query spectrum to existing library
#'
#' @export
#'
#' @importFrom MSnbase fData readMgfData
#' @importFrom tools file_ext
#' @importFrom stringr str_replace_all fixed
#' @importFrom OrgMassSpecR SpectrumSimilarity

library_similarity<-function(library, query_spectrum=NULL , method = c("Simple", "Cosinus"), tops = 3, prec_mz = 1000, ppm_search = 20, relative = 1){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  #################
  ### Check inputs:
  #################

  if (missing(library)){
    stop("Please provide the output of library_generator() or a .mgf file as input library!")}

  if (is.character(library)){
    if (file_ext(library)!="mgf"){
      stop("The file extension of your input library must be mgf!")
    }}

  if (is.list(library)){
    if (length(library)==2 & "complete" %in% names(library2)){
      library = library$complete
    }
    if (length(library)!=2 || (!is.list(library$sp)) || !is.data.frame(library$metadata)){
      stop("Please make sure your input library is a valid output of library_generator()!")
    }}

  if (is.null(query_spectrum)){
      stop("Please provide a 2 column query spectrum!")
  } else {
  if (ncol(query_spectrum)<2){
    stop("Spectrum must have 2 columns m/z and intensity!")
    }
  }

  ###############################
  ### Preprocess query spectrum:
  ###############################

  dat = query_spectrum[,1:2]

  # Normalize, cut only masses smaller than precursor and filter background noise:

  dat[,2]=dat[,2]/max(dat[,2])*100
  selected = which((dat[,1] < prec_mz+10) & (dat[,2]>relative))

  if (length(selected)>0){
      dat = dat[selected,]
    dat = data.matrix(matrix(dat,ncol=2))
    abs_search = ppm_search/1000000*median(dat[,1]) # Da window for spectra search
  } else {
    stop("Spectrum not valid!")
  }

  #####################################
  ### Reading from spectral library:
  #####################################

  if (is.character(library)){ # If input is a mgf file name
    library=readMGF2(library)}

  metadata = library$metadata
  spectrum_list = library$sp

  ###########################
  ### Run spectra search:
  ###########################

  NS = nrow(metadata)
  tops = min(tops, NS)

  nb_matched_scores = rep(0,NS)
  sim_scores = rep(0,NS)

  for (i in 1:NS){

    if (method == "Simple"){
      dist_spec = abs(sapply(spectrum_list[[i]][,1], "-", dat[,1]))
      sim_scores[i] = sum(dist_spec<=abs_search)
     }

    if (method == "Cosinus"){
      sim_scores[i] = SpectrumSimilarity(dat, spectrum_list[[i]], t = abs_search)
    }
    graphics.off()
  }

  ###########################
  ### Plot and output results:
  ###########################

  # Top scores and output:

  indexes = order(sim_scores,decreasing=T)[1:tops]
 # indexes = intersect(indexes,which(sim_scores>0))

  SELECTED_LIBRARY = list()
  SELECTED_LIBRARY$sp = library$sp[indexes]
  SELECTED_LIBRARY$metadata = library$metadata[indexes,]

  print(SELECTED_LIBRARY$sp)
  print(indexes)
  # Plot:
  for (i in 1:length(indexes)){
    bottom.label = paste0("Library spectrum of ID = ", SELECTED_LIBRARY$metadata$ID[i])
    xlim = c(min(dat[,1])*0.8, max(dat[,1])*1.2)
    SpectrumSimilarity(dat,  SELECTED_LIBRARY$sp[[i]], t = abs_search,
                       top.label="Query spectrum",bottom.label= bottom.label,xlim=xlim)
  }
  return(SELECTED_LIBRARY)
}
