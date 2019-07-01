#' Create pseudo-metadata for library generation
#'
#' Function used by library_generator to extract LC-MS features from MS1 chromatogram of DDA/targeted-mode
#'
#' @importFrom xcms xcmsSet
#' @importFrom CAMERA xsAnnotate groupFWHM findIsotopesWithValidation getPeaklist
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
metadata_MS1_screener<-function(raw_data_file, ref = NULL, max.charge = 1,
                                ppm_search = 20, baseline = 1000, snthreshold = 30, min.points = 8){

  options(stringsAsFactors = F)
  options(warn=-1)

  ####################################
  ### Extract features from MS1 scans
  ####################################

  xs <- xcmsSet(raw_data_file, method = "centWave", ppm = ppm_search*2,
                snthresh = snthreshold, prefilter = c(min.points, baseline))

  ### Remove isotopes

  an <- xsAnnotate(xs)
  an <- groupFWHM(an, perfwhm=2)
  an <- findIsotopesWithValidation(an, ppm = ppm_search,  mzabs = 0.01, maxcharge = max.charge, database = "pubchem")  # optional but recommended.

  peaks = getPeaklist(an)
  peaks = data.frame(peaks)

  iso_pos = grep("[M]+", peaks$isotopes, fixed=T)
  iso_neg = grep("[M]-", peaks$isotopes, fixed=T)

  if (length(iso_pos)>0){
    polarity = "Positive"
    adduct = "M+H"
    iso1 = grep("[M+1]+", peaks$isotopes, fixed = T)
    iso2 = grep("[M+2]+", peaks$isotopes, fixed = T)
    iso3 = grep("[M+3]+", peaks$isotopes, fixed = T)
    iso4 = grep("[M+4]+", peaks$isotopes, fixed = T)
  } else {
    polarity = "Negative"
    adduct = "M-H"
    iso1 = grep("[M+1]-", peaks$isotopes, fixed=T)
    iso2 = grep("[M+2]-", peaks$isotopes, fixed= T)
    iso3 = grep("[M+3]-", peaks$isotopes, fixed= T)
    iso4 = grep("[M+4]-", peaks$isotopes, fixed= T)
  }

  peaks = peaks[-c(iso1, iso2, iso3, iso4),]

  ###########################################
  ### Transform peak table to usable metadata
  ###########################################

  ID_list = paste0(file_path_sans_ext(raw_data_file),"_",1:nrow(peaks))
  MS1_metadata = cbind(peaks$mz, peaks$rt/60, polarity, adduct, 1, ID_list)
  MS1_metadata = data.frame(MS1_metadata)
  colnames(MS1_metadata) = c("PEPMASS","RT","IONMODE","ADDUCT","CHARGE","ID")
  MS1_masses = peaks$mz

  ##################################
  ### Filter user provided metadata
  #################################

  if (!is.null(ref)){
    ref1 = c() # New metadata
    for (i in 1:nrow(ref)){
      ppm_errors = ppm_distance(MS1_masses,ref$PEPMASS[i])
      valid = which(ppm_errors<=ppm_search & MS1_metadata$ADDUCT == ref$ADDUCT[i])
      nbf = length(valid)
      if (nbf>0){ # Update the found feature:
        ref_selected = do.call("rbind", replicate(nbf, ref[i,], simplify = FALSE))
        ref_selected$RT = MS1_metadata$RT[valid]
        if (nbf>1){ref_selected$ID = paste0(ref_selected$ID,"_",1:nbf)}
        ref1 = rbind.data.frame(ref1,ref_selected)
      }
    }
  } else {
    ref1 = MS1_metadata
  }
  return(ref1)
}

############################
### Internal functions:
###########################

ppm_distance<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  return(abs((x-y)/y*1000000))
}
