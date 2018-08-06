#' Read and combine targeted MS1 scans from one LC-MS/MS file
#'
#' Function used by library_generator to detect MS1 scans
#' @export

process_MS1<-function(mzdatafiles,ref,rt_search=10,ppm_search=20,MS2_type = c("DDA","Targeted"),baseline= 1000,normalized=T){

  ### Initialize variables

  MS1_Janssen = NULL
  new_MS1_meta_data = c() # Compounds detected in the ms data
  MS1_scan_list = list() # List of spectrum2 objects
  scan_number = c() # Which scan number in the raw chromatogram
  new_PEP_mass = c() # The real mass in samples
  mass_dev = c() # Mass deviations in ppm
  spectrum_list = list() # List of spectra to save
  N=0

  ### Read the raw data file

  MS1_Janssen <- try(readMSData(mzdatafiles, msLevel = 1, verbose = FALSE),silent=T)

  if (class(MS1_Janssen)!="try-error"){ # If data contains MS1 scan

    print(paste0("Processing MS1 scans of data file ",mzdatafiles," ..."))

    ### Extract useful information:

    prec_theo=ref$PEPMASS
    prec_rt=as.numeric(ref$RT)*60 # Allow N/A
    MS1_prec_rt = rtime(MS1_Janssen)
    int_max_list = c() # Maximal intensity of MS1 mass

    for (i in 1:nrow(ref)){

    # Save highest intensity MS1 data:

      if (!is.na(prec_rt[i])){
          scan_range=which(MS1_prec_rt >= prec_rt[i] - rt_search & MS1_prec_rt <= prec_rt[i] + rt_search)
     } else {scan_range=1:length(MS1_prec_rt)}

    # Find the scan that corresponds to the meta data

    int_max=0
    valid_k=0

    if (length(scan_range)>0){
      for (k in scan_range){ # Check whether the precursor peak is detected
        Frag_data = MS1_Janssen[[k]]

        if (length(Frag_data@mz)>1){ # At least the scan is not empty

        error = min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
        valid = which.min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
        int = Frag_data@intensity[valid] # Precursor intensity

      # We now check whether the precursor is precisely isolated
        if (error<=ppm_search & int>int_max){
          valid_k=k
          int_max=int}
      }}}

    if ((valid_k!=0) & (int_max > baseline*5)) { # If the scan is found and signal higher than 5 times of baseline
        scan_number = c(scan_number,valid_k)  # Save scan number
        int_max_list=c(int_max_list, int_max) # Save maximal intensity

        masslist=MS1_Janssen[[valid_k]]@mz
        wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
        dev_ppm= min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
        dev_ppm=round(dev_ppm,2)

        new_PEP_mass = c(new_PEP_mass,masslist[wp]) # Save detected mass
        mass_dev = c(mass_dev,dev_ppm)

        MS1_scan_list[[N+1]]=MS1_Janssen[[valid_k]]
        N=N+1
        new_MS1_meta_data = rbind(new_MS1_meta_data,ref[i,])
    }}

  ### Update metadata

    new_MS1_meta_data[,"PEPMASS"]=new_PEP_mass
    new_MS1_meta_data[,"RT"]= round(MS1_prec_rt[scan_number]/60,2)
    new_MS1_meta_data[,"FILENAME"]=rep(mzdatafiles,N)
    new_MS1_meta_data[,"MSLEVEL"]=rep(1,N)
    new_MS1_meta_data[,"TIC"]= int_max_list
    new_MS1_meta_data[,"PEPMASS_DEV"]=mass_dev
    new_MS1_meta_data[,"SCAN_NUMBER"] = scan_number

  ### Denoise spectra

  if (!is.null(new_MS1_meta_data)){
    included=c()
    n0=0

    for (i in 1:N){
      dat = cbind(MS1_scan_list[[i]]@mz,MS1_scan_list[[i]]@intensity)
      # Filter background noise and choose only peaks until +7Da of the precursor!!
      selected = which(dat[,1]>new_MS1_meta_data$PEPMASS[i]-0.5 & dat[,1] < new_MS1_meta_data$PEPMASS[i]+7 & dat[,2]>baseline)
      if (length(selected)>2){ # At least 3 peaks should be present for isotope determination
        dat = dat[selected,]
        dat = matrix(dat,ncol=2)
        if (normalized){dat[,2]=dat[,2]/max(dat[,2])*100}
        n0=n0+1
        spectrum_list[[n0]]=dat
        included= c(included,i)}}

    ### Keep only no empty spectra
    new_MS1_meta_data = new_MS1_meta_data[included,]
  }
  } else {
   print(paste0("No MS1 scan in the data file ",mzdatafiles," !"))
 }

  return(list(sp=spectrum_list,metadata=new_MS1_meta_data))
}

