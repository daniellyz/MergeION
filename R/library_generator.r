#' Batch processing of LC-MS/MS chromatograms into a spectral library
#'
#' The function picks up scans according to m/z (and retention time) specified in the metadata and merge them into a spcetral library (new or existing). The raw LC-MS/MS files must be centroid-mode mzML, mzMXL or mzData
#'
#' @param raw_data_files A character vector of file names of chromatograms from which scans are extracted. All files must have be in centroid-mode with mzML or mzMXL extension!
#' @param metadata_file File name of the metadata. Must be a single character with csv extension. The first five columns of the metadata must be (in order): "PEPMASS" (precursor masses that we want to find in chromatograms), "RT" (retention time of metabolic features to be found, in minute, please put it to N/A if unknown), "IONMODE" (must be "Positive" or "Negative"),"ADDUCT" (precursor ion adduct type, must be one of "M+H","M+Na","M+K","M-H" and "M+Cl") and "ID" (A unique identifier for targeted compounds in spectral library).
#' @param mslevel Must be 1 (if only MS1 scans/isotopic patterns of targeted m/z are extracted), 2 (if only MS2 scans are extracted) or c(1,2) (if both MS1 and MS2 scans are extracted). Note: Isotopic patterns in MS1 scans are useful for determining precursor formula !
#' @param MS2_type  A single character ("DDA" or "Targeted") if all raw_dat_files are acquired in the same mode; A character vector precising the acquisition mode of each file in raw_data_files (e.g. c("DDA","Targeted","DDA"))
#' @param rt_search Retention time search tolerance (in second) for targeted RT
#' @param ppm_search m/z search tolerance (in ppm) for targeted m/z
#' @param baseline Numeric if all raw_dat_files have the same beseline intensity (the minimum intensity that is considered as a mass peak and written into the library); or a numeric vector describing the baseline of each file (e.g. c(100,4000,10))
#' @param normalized Logical, TRUE if the intensities of extracted spectra need to normalized so that the intensity of highest peak will be 100
#' @param input_library Name of the library into which new scans are added, the file extension must be mgf; please set to empty string "" if the new library has no dependency with previous ones.
#' @param output_library Name of the output library, the file extension must be mgf
#'
#' @return
#' \itemize{
#'   \item{"library$sp" ~ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{"library$metadata" ~ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following five columns are added: FILENAME, MSLEVEL, TIC, MASS_DEV, SCANNUMBER and SCANS}
#'   \item{"<ouput_library>" ~ A mgf spectral library file will be written in user's working directory. It contains both spectra and metadata}
#'   \item{"<ouput_library.txt>" ~ Metadata will be written as a tab-seperated .txt file in user's working directory. Users can check this file in excel or open office.}
#' }
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples

#' ### We download our three test data sets:
#' # Details of these data can be found at: https://zenodo.org/record/1322562
#'
#' url = "https://zenodo.org/record/1326555/files/"
#' original_files = c("TESTMIX2_180504_MAS011_06.mzXML",
#'                    "JNJ42165279_171214_MAS006_14.mzXML",
#'                    "GMP_R601592_150925_MAS006_04.mzXML")
#' download.file(paste0(url,original_files[1]),destfile="F1.mzXML") # Download and rename the files
#' download.file(paste0(url,original_files[2]),destfile="F2.mzXML")
#' download.file(paste0(url,original_files[3]),destfile="F3.mzXML")
#'
#' ### It is time to batch-process the first two files and create our first spectral library:
#'
#' raw_data_files = c("F1.mzXML","F2.mzXML")
#' metadata_file = paste0(url,"library_metadata.csv")
#' mslevel = c(1,2)
#' # Both MS1 and MS2 scans are extracted! MS1 scan contains isotopic pattern of targeted m/z
#' MS2_type = c("DDA","Targeted") # Mode of MS/MS experiment for F1 and F2 respectively
#' rt_search = 12 # Retention time tolerance (s)
#' ppm_search = 10  # Mass tolerance (ppm)
#' baseline = 1000 # Noise level of each mass spectrum, 1000 is fixed for both chromatograms
#' input_library = "" # A brand new library, there's no previous dependency
#' output_library = "library_V1.mgf" # Name of the library
#'
#' library1 = library_generator(raw_data_files, metadata_file, mslevel, MS2_type, rt_search, ppm_search,
#'       baseline,normalized = T, input_library, output_library)
#'
#' ### We added the targeted scans of F3.mzXML to spectral library version 2:
#'
#' raw_data_files = "F3.mzXML" # The new LC-MS/MS data
#' MS2_type = "DDA"
#' input_library = "library_V1.mgf" # The first mgf file of library1
#' output_library = "library_V2.mgf" # The name of the new spectral library
#'
#' library2 = library_generator(raw_data_files, metadata_file, mslevel, MS2_type, rt_search, ppm_search,
#'       baseline, normalized = T, input_library, output_library)
#'
#' ### In the end, two spectral library versions "library_V1.mgf" and "library_V2.mgf" should appear in the working directory along with metadata table (txt files)
#'
#' @export
#'
#' @importFrom MSnbase readMSData rtime tic fData readMgfData precursorMz
#' @importFrom tools file_ext
#' @importFrom utils write.table read.csv

library_generator<-function(raw_data_files,metadata_file,mslevel = c(1,2),MS2_type = "DDA",rt_search = 12,ppm_search = 20,
                            baseline = 1000, normalized=T,
                            input_library="", output_library=""){

  options(stringsAsFactors = FALSE)
  options(warn=-1)
  FF = length(raw_data_files)
  spectrum_list = list()
  metadata = c()
  NN = 0
  unlink(output_library)

  ##############################
  ### Check function inputs:
  ##############################

  if (missing(raw_data_files) || (!is.vector(raw_data_files))){
    stop("Please provide a list of chromatogram files!")}

  if (!all(file_ext(raw_data_files) %in% c("mzML","mzXML","mzData"))){
    stop("Chromatogram files must be in mzML, mzXML or mzData format!")}

  if (missing(metadata_file) || (!is.character(metadata_file)) || (length(metadata_file)!=1) || file_ext(metadata_file)!="csv"){
    stop("Metadata must be written in one single csv!")}

  if (!all(mslevel %in% c(1,2))){
    stop("mslevel must be 1 or 2!")}

  if (length(MS2_type)==1){
    MS2_type = rep(MS2_type,FF)
  } else {
    if (length(MS2_type)!=FF){
      stop("The length of MS2_type must be the same as raw_data_files!")}}

  if (length(baseline)==1){
    baseline = rep(baseline,FF)
  } else {
    if (length(baseline)!=FF){
      stop("The length of baseline must be the same as raw_data_files!")}}

  MS2_type = match.arg(MS2_type,choices=c("DDA","Targeted"),several.ok = TRUE)

  if (input_library!=""){
    if (file_ext(input_library)!="mgf"){
      stop("The input library must be mgf format!")
    }}

  if (missing(output_library) || file_ext(output_library)!="mgf"){
    stop("The output library in mgf format must be filled!")
  }

  if (input_library==output_library){
    stop("The new library must be saved under a different name as the previous library!")
  }

  #######################################
  ### Read from metadata and old library:
  #######################################

  ref<-read.csv(metadata_file,sep=";",dec=".",header=T,stringsAsFactors = F)
  if (is.null(nrow(ref))){
    labels=names(ref)
    ref=data.frame(matrix(ref,nrow=1))
    colnames(ref)=labels}

  ref[,5]=as.character(ref[,5]) # Make sure IDs are characters
  colnames(ref)[1]="PEPMASS"  # Make sure column name is correct!
  colnames(ref)[2]="RT"
  colnames(ref)[3]="IONMODE"
  colnames(ref)[4]="ADDUCT"
  colnames(ref)[5]="ID"

  for (j in 1:ncol(ref)){  # Change column names to avoid duplicates
    if (colnames(ref)[j] %in% c("FILENAME","MSLEVEL","TIC","PEPMASS_DEV","PEPMASS_DEV","SCAN_NUMBER","SCANS")){
      colnames(ref)[j]=paste0(colnames(ref)[j],"_000")
    }}

  if (ncol(ref)<5){
    stop("Metadata must contain at least 5 columns!")}

  if (!is.numeric(ref$PEPMASS)){
    stop("Precursor masses (PEPMASS) must be numeric!")}

  if (!all(ref$IONMODE %in% c("Positive","Negative"))){
    stop("Ion mode must be Postive or Negative!")}

  if (!all(ref$ADDUCT %in% c("M+H","M+Na","M+K","M-H","M+Cl"))){
    stop("Adduct type not valid!")}

  if (input_library!=""){
    old_dat=readMgfData(input_library, verbose = FALSE)
    spectrum_list=Mgf2Splist(old_dat)
    metadata=fData(old_dat)
    metadata_items=colnames(metadata)[1:(ncol(metadata)-6)]
    if (!(all(colnames(ref)==metadata_items))){
      stop("The new library must contain same exact metadata as the previous one!")}
    metadata=metadata[,1:(ncol(metadata)-1)] # Remove SCANS!
    NN = length(spectrum_list)
  }

  ##################
  ### Batch process:
  ##################

  for (ff in 1:FF){

    if (2 %in% mslevel){
      dat2 = process_MS2(raw_data_files[ff],ref,rt_search,ppm_search, MS2_type[ff], baseline[ff])
      LL2 = length(dat2$sp) # Added library size
      ref2 = dat2$ref_MS2 # Filter metadata data for MS1 searcch
      if (LL2>0){
        for (n in 1:LL2){spectrum_list[[NN+n]]=dat2$sp[[n]]} # Update spectrum list
        metadata=rbind(metadata,dat2$metadata) # Update metadata
        NN=NN+LL2
      }}

    if (1 %in% mslevel){
      dat1 = process_MS1(raw_data_files[ff],ref2,rt_search,ppm_search, baseline[ff])
      LL1= length(dat1$sp) # Added library size
      if (LL1>0){
        for (n in 1:LL1){spectrum_list[[NN+n]]=dat1$sp[[n]]} # Update spectrum list
        metadata=rbind(metadata,dat1$metadata) # Update metadata
        NN=NN+LL1
      }}
  }

  ####################
  ### Return results:
  ####################

  library = list()
  library$sp = spectrum_list
  library$metadata = cbind.data.frame(metadata,SCANS=1:nrow(metadata))

  writeMGF2(library$sp,library$metadata,output_library)
  write.table(library$metadata,paste0(output_library,".txt"),col.names = T,row.names=F,dec=".",sep="\t")
  return(library)
}

###########################
### Internal functions:
###########################

Mgf2Splist<-function(MGFdat){

  # From a MSnBase object to a list of spectra m/z intensity
  N=length(MGFdat)
  spectrum_list=list()
  for (i in 1:N){spectrum_list[[i]]=cbind(MGFdat[[i]]@mz,MGFdat[[i]]@intensity)}
  return(spectrum_list)
}

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
