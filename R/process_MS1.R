process_MS1<-function(mzdatafiles,ref,rt_search=10,ppm_search=20,search_adduct = T,baseline= 1000,normalized=T){

  # The function reads from 1 raw data and extracts a list of targeted spectra (Scans) and modified metadata

  ### Read new spectra and meta data:

  prec_theo=ref$PEPMASS
  prec_rt=as.numeric(ref$RT)*60 # Allow N/A

  ### Initialize variables

  MS1_Janssen = NULL
  new_MS1_meta_data = c() # Compounds detected in the ms data
  MS1_scan_list = list() # List of spectrum2 objects
  scan_number = c() # Which scan number in the raw chromatogram is found
  new_PEP_mass = c() # The real mass in samples
  ADDUCT_LABEL = c()
  spectrum_list = list() # List of spectra to save
  N=0

  ### Read the raw data file

  MS1_Janssen <- readMSData(mzdatafiles, msLevel = 1, verbose = FALSE)

  if (length(MS1_Janssen)>0){ # If data contains MS1 scan

    MS1_prec_rt = rtime(MS1_Janssen)
    int_max_list = c() # Maximal intensity of MS1 mass

    for (i in 1:nrow(ref)){

    # Save highest intensity MS1 data:

      if (!is.na(prec_rt[i])){
          scan_range=which(MS1_prec_rt >= prec_rt[i] - rt_search & MS1_prec_rt <= prec_rt[i] + rt_search)
     } else {scan_range=1:length(MS1_prec_rt)}

    # Calculate adducts

    if (search_adduct & (ref$IONMODE[i]=="Positive")){
      prec_adduct = prec_theo[i] - 1.007276 + 22.989221} # Soidum

    if (search_adduct & (ref$IONMODE[i]=="Negative")){
      prec_adduct = prec_theo[i] + 1.007276 + 34.969401} # Chlore

    # Find the scan that corresponds to the meta data

    int_max=0
    valid_k=0
    new_adduct_type=0

    if (length(scan_range)>0){
      for (k in scan_range){ # Check whether the precursor peak is detected
        Frag_data = MS1_Janssen[[k]]
        adduct_type = 0 # Default adduct type

        error= min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
        valid= which.min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
        if (search_adduct){
          error1 = min(abs(Frag_data@mz-prec_adduct)/prec_adduct*1000000)
          valid1 = which.min(abs(Frag_data@mz-prec_adduct)/prec_adduct*1000000)
          if (error1<error){ # If adduct error smaller than original
            error = error1
            valid = valid1
            adduct_type = 1}}

      # We now check whether the precursor is precisely isolated
        if (error<=ppm_search & Frag_data@intensity[valid]>int_max){
          valid_k=k
          int_max=Frag_data@intensity[valid]
          new_adduct_type=adduct_type} # To check whether is original or adduct detected
      }}

    if ((valid_k!=0) & (int_max > baseline*5)) { # If the scan is found and signal higher than 5 times of baseline
        scan_number = c(scan_number,valid_k)  # Save scan number
        int_max_list=c(int_max_list, int_max) # Save maximal intensity

        masslist=MS1_Janssen[[valid_k]]@mz
        if (new_adduct_type==0){wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)}
        if (new_adduct_type==1){wp= which.min(abs(masslist-prec_adduct)/prec_adduct*1000000)}
        new_PEP_mass = c(new_PEP_mass,masslist[wp]) # Save detected mass

        if (ref$IONMODE[i]=="Positive" & new_adduct_type==1){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+Na")}
        if (ref$IONMODE[i]=="Positive" & new_adduct_type==0){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+H")}
        if (ref$IONMODE[i]=="Negative" & new_adduct_type==1){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+Cl")}
        if (ref$IONMODE[i]=="Nagative" & new_adduct_type==0){ADDUCT_LABEL=c(ADDUCT_LABEL,"M-H")}

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
    new_MS1_meta_data[,"ADDUCT"]=ADDUCT_LABEL
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
  }}
  return(list(sp=spectrum_list,metadata=new_MS1_meta_data))
}

