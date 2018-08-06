#' Read and combine targeted MS2 scans from one LC-MS/MS file
#'
#' Function used by library_generator to detect MS2 scans
#' @export

process_MS2<-function(mzdatafiles,ref,rt_search=10,ppm_search=20, MS2_type = c("DDA","Targeted"),baseline= 1000,normalized=T){

  ### Initialize variables

  MS2_Janssen = NULL
  new_MS2_meta_data = c() # Compounds detected in the ms data
  MS2_scan_list = list() # List of spectrum2 objects
  scan_number = c() # Which scan number in the raw chromatogram is found
  new_PEP_mass = c() # The real mass in samples
  mass_dev = c() # Mass deviations in ppm
  spectrum_list = list() # List of spectra to save
  N=0 # Numerator

  ### Read the raw data file

  MS2_Janssen <- try(readMSData(mzdatafiles, msLevel = 2,  verbose = FALSE),silent=T)

  if (class(MS2_Janssen)!="try-error"){ # If data contains MS2 scan

    print(paste0("Processing MS2 scans of data file ",mzdatafiles," ..."))

    ### Filter ref because not all targeted m/z exists or fragmented in the sample! Important for not to search whatever

    targets =  unique(precursorMz(MS2_Janssen))

    dev_targets = sapply(ref$PEPMASS,function(x) min(abs(x-targets)))
    valid = which(dev_targets <= 1) #  Find targeted metadata in experimental file!!
    ref = ref[valid,]

    if (nrow(ref)>0){

    ### Extract useful informations

      prec_theo=ref$PEPMASS
      prec_rt=as.numeric(ref$RT)*60 # Allow N/A
      MS2_prec_rt = rtime(MS2_Janssen) # In second
      MS2_tic = tic(MS2_Janssen)

    ### Check one by one targeted m/z:

      for (i in 1:nrow(ref)){

      # Search RT window if provided:

        if (!is.na(prec_rt[i])){
        scan_range=which(MS2_prec_rt >= prec_rt[i] - rt_search & MS2_prec_rt <= prec_rt[i] + rt_search)
        } else {
        scan_range=1:length(MS2_prec_rt)}

     # Find the scan that corresponds to the meta data

      tic_max=0
      valid_k=0

      if (length(scan_range)>0){
        # Check scan by scan the masses:
        for (k in scan_range){

          Frag_data = MS2_Janssen[[k]]

          if (length(Frag_data@mz)>1){ # At least the scan is not empty

             precursor_abundant = 1 # Whether the "precursor" found is abundant

             if (MS2_type=="DDA"){
                 error=abs(Frag_data@precursorMz-prec_theo[i])/prec_theo[i]*1000000}

             if (MS2_type=="Targeted"){ # Some issues with certain instrument that targeted scans are badly labeled (not the real precursor)
                 error= min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
                 prec_ind = which.min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)

                 # Make sure that the precursor mass detected in targeted mode must be at least 5 times higher than baseline!
                 if (Frag_data@intensity[prec_ind]<baseline*5) {
                    precursor_abundant=0} else {precursor_abundant=1}
             }

      # We now check whether the isolated scan is better than previous:
            if ((error<=ppm_search) & MS2_tic[k]>tic_max & precursor_abundant==1){
              valid_k=k # Update the scan number
              tic_max=MS2_tic[valid_k]}
      }}}

    if (valid_k==0){}

    if (valid_k!=0){ # If the scan is found

      scan_number = c(scan_number,valid_k)  # Save scan number

      ### Update metadata:
      # Save detected precursor mass:

        if (MS2_type=="DDA"){
          mz = MS2_Janssen[[valid_k]]@precursorMz
          dev_ppm= min(abs(mz-prec_theo[i])/prec_theo[i]*1000000)
          dev_ppm=round(dev_ppm,2)
          new_PEP_mass = c(new_PEP_mass,mz)
          mass_dev = c(mass_dev,dev_ppm)
        }

        if (MS2_type=="Targeted"){
          masslist=MS2_Janssen[[valid_k]]@mz
          wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
          dev_ppm= min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
          dev_ppm=round(dev_ppm,2)
          new_PEP_mass = c(new_PEP_mass,masslist[wp])
          mass_dev = c(mass_dev,dev_ppm)
        }

        MS2_scan_list[[N+1]]=MS2_Janssen[[valid_k]]
        N=N+1
        new_MS2_meta_data = rbind(new_MS2_meta_data,ref[i,])
    }}

  ### Update metadata
    new_MS2_meta_data[,"PEPMASS"]=new_PEP_mass
    new_MS2_meta_data[,"RT"]= round(MS2_prec_rt[scan_number]/60,2)  # minutes
    new_MS2_meta_data[,"FILENAME"]=rep(mzdatafiles,N)
    new_MS2_meta_data[,"MSLEVEL"]=rep(2,N)
    new_MS2_meta_data[,"TIC"]= MS2_tic[scan_number]
    new_MS2_meta_data[,"PEPMASS_DEV"]=mass_dev
    new_MS2_meta_data[,"SCAN_NUMBER"] = scan_number

  ### Denoise spectra
  if (!is.null(new_MS2_meta_data)){
    included = c() # Not filtered
    n0=0

    for (i in 1:N){
      dat = cbind(MS2_scan_list[[i]]@mz,MS2_scan_list[[i]]@intensity)

      # Cut only masses smaller than precursor and filter background noise:
      selected = which((dat[,1] < new_MS2_meta_data$PEPMASS[i]+10) & (dat[,2]>baseline))
      if (length(selected)>0){
        dat = dat[selected,]
        dat = matrix(dat,ncol=2)
        if (normalized){dat[,2]=dat[,2]/max(dat[,2])*100}
        n0=n0+1
        spectrum_list[[n0]]=dat
        included= c(included,i)}}

  ### Keep only no empty spectra
      new_MS2_meta_data = new_MS2_meta_data[included,]
      id_kept = new_MS2_meta_data$ID
      ref = ref[match(id_kept,ref$ID),]
     }
    } else {
      print(paste0("No MS2 scan in the data file ",mzdatafiles," matches with metadata!"))}
   } else {
  print(paste0("No MS2 scan in the data file ",mzdatafiles," !"))
 }

  return(list(sp=spectrum_list,metadata=new_MS2_meta_data,ref_MS2=ref))
}

