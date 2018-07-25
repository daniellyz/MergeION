#' Read and combine targeted MS2 scans from one LC-MS/MS file
#'
#' Function used by library_generator to detect MS2 scans

process_MS2<-function(mzdatafiles,ref,rt_search=10,ppm_search=20,search_adduct = T,MS2_type = c("DDA","Targeted"),baseline= 1000,normalized=T){


  ### Read new spectra and meta data:

  prec_theo=ref$PEPMASS
  prec_rt=as.numeric(ref$RT)*60 # Allow N/A

  ### Initialize variables

  MS2_Janssen = NULL
  new_MS2_meta_data = c() # Compounds detected in the ms data
  MS2_scan_list = list() # List of spectrum2 objects
  scan_number = c() # Which scan number in the raw chromatogram is found
  new_PEP_mass = c() # The real mass in samples
  ADDUCT_LABEL = c() # Type of adduct of each scan
  spectrum_list = list() # List of spectra to save
  N=0 # Numerator

  ### Read the raw data file

  MS2_Janssen <- readMSData(mzdatafiles, msLevel = 2,  verbose = FALSE)

  if (length(MS2_Janssen)>0){ # If data contains MS2 scan

    MS2_prec_rt = rtime(MS2_Janssen) # In second
    MS2_tic = tic(MS2_Janssen)

    ### Check one by one targeted m/z:

    for (i in 1:nrow(ref)){

      # Search RT window if provided:

      if (!is.na(prec_rt[i])){
        scan_range=which(MS2_prec_rt >= prec_rt[i] - rt_search & MS2_prec_rt <= prec_rt[i] + rt_search)
      } else {
        scan_range=1:length(MS2_prec_rt)}

     # Calculate adducts

     if (search_adduct & (ref$IONMODE[i]=="Positive")){
       prec_adduct = prec_theo[i] - 1.007276 + 22.989221} # Soidum

     if (search_adduct & (ref$IONMODE[i]=="Negative")){
       prec_adduct = prec_theo[i] + 1.007276 + 34.969401} # Chlore

     # Find the scan that corresponds to the meta data

      tic_max=0
      valid_k=0
      new_adduct_type=0

      if (length(scan_range)>0){
        # Check scan by scan the masses:
        for (k in scan_range){
          Frag_data = MS2_Janssen[[k]]
          adduct_type = 0 # Default is no adduct!

          if (MS2_type=="DDA"){
            error=abs(Frag_data@precursorMz-prec_theo[i])/prec_theo[i]*1000000
            if (search_adduct){
              error1 = abs(Frag_data@precursorMz-prec_adduct)/prec_adduct*1000000
              if (error1<error){ # If adduct error smaller than original
                error = error1
                adduct_type = 1}}}

          if (MS2_type=="Targeted"){
            error= min(abs(Frag_data@mz-prec_theo[i])/prec_theo[i]*1000000)
          if (search_adduct){
            error1 = min(abs(Frag_data@mz-prec_adduct)/prec_adduct*1000000)
            if (error1<error){ # If adduct error smaller than original
              error = error1
              adduct_type = 1}}}

      # We now check whether the isolated scan is better than previous:
          if ((error<=ppm_search) & MS2_tic[k]>tic_max){
            valid_k=k
            tic_max=MS2_tic[valid_k]
            new_adduct_type=adduct_type} # To check whether is original or adduct detected
      }}

    if (valid_k==0){}

    if (valid_k!=0){ # If the scan is found

      scan_number = c(scan_number,valid_k)  # Save scan number

      ### Update metadata:

        if (MS2_type=="DDA"){new_PEP_mass = c(new_PEP_mass,MS2_Janssen[[valid_k]]@precursorMz)}
        if (MS2_type=="Targeted"){
          masslist=MS2_Janssen[[valid_k]]@mz
          if (new_adduct_type==0){wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)}
          if (new_adduct_type==1){wp= which.min(abs(masslist-prec_adduct)/prec_adduct*1000000)}
          new_PEP_mass = c(new_PEP_mass,masslist[wp])} # Save detected precursor mass

        if (ref$IONMODE[i]=="Positive" & new_adduct_type==1){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+Na")}
        if (ref$IONMODE[i]=="Positive" & new_adduct_type==0){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+H")}
        if (ref$IONMODE[i]=="Negative" & new_adduct_type==1){ADDUCT_LABEL=c(ADDUCT_LABEL,"M+Cl")}
        if (ref$IONMODE[i]=="Nagative" & new_adduct_type==0){ADDUCT_LABEL=c(ADDUCT_LABEL,"M-H")}

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
    new_MS2_meta_data[,"ADDUCT"] = ADDUCT_LABEL
    new_MS2_meta_data[,"SCAN_NUMBER"] = scan_number

  ### Denoise spectra
  if (!is.null(new_MS2_meta_data)){
    included = c() # Not filtered
    n0=0

    for (i in 1:N){
      dat = cbind(MS2_scan_list[[i]]@mz,MS2_scan_list[[i]]@intensity)

      # Cut only masses smaller than precursor and filter background noise:
      selected = which((dat[,1] < new_MS2_meta_data$PEPMASS[i]+3) & (dat[,2]>baseline))
      if (length(selected)>0){
        dat = dat[selected,]
        dat = matrix(dat,ncol=2)
        if (normalized){dat[,2]=dat[,2]/max(dat[,2])*100}
        n0=n0+1
        spectrum_list[[n0]]=dat
        included= c(included,i)}}

  ### Keep only no empty spectra
    new_MS2_meta_data = new_MS2_meta_data[included,]
  }}
  return(list(sp=spectrum_list,metadata=new_MS2_meta_data))
}
