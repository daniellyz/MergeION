library(MSnbase)
library(tools)
source("process_MS1.R")
source("process_MS2.R")
source("writeMGF2.R")

options(stringsAsFactors = FALSE)
options(warn=-1)

library_writer_MS1<-function(raw_data_files,metadafile,output_library,rt_search = 12,ppm_search = 20,search_adduct = T,baseline = 1000, normalized=T, old_library=NULL){
  
  # This function batch processes raw data files, and write them into a MGF file (MS1). It also takes into account historical library
  
  spectrum_list=list()
  metadata=c()
  NN=0 
  
  # Read from previous library:
  
  if (!is.null(old_library)){
    old_dat=readMgfData(old_library, verbose = FALSE)
    spectrum_list=Mgf2Splist(old_dat) 
    metadata=fData(old_dat)
    metadata=metadata[,1:(ncol(metadata)-1)]
    NN = length(spectrum_list)
  }
  
  # Batch process:
  
  FF = length(raw_data_files)
  
  for (ff in 1:FF){
    
    dat1 = process_MS1(raw_data_files[ff],metadafile,rt_search,ppm_search,search_adduct, baseline)
    nls= length(dat1$sp) # Added library size
      if (length(nls)>0){
      for (n in 1:nls){spectrum_list[[NN+n]]=dat1$sp[[n]]} # Update spectrum list
      metadata=rbind(metadata,dat1$metadata) # Update meta data
      NN = NN+nls} # Update size
  }
  
  # Filter: only keep the highest TIC
  valid=filter_scans(metadata)
  metadata = metadata[valid,]
  metadata = cbind.data.frame(metadata,SCANS=1:nrow(metadata))
  spectrum_list = spectrum_list[valid]
  
  # Save results:
  write.table(metadata,paste0(output_library,".txt"),col.names = T,row.names=F,dec=".",sep="\t")
  save(spectrum_list,metadata,file=paste0(output_library,".rData"))
  writeMGF2(spectrum_list,metadata,output_library)
}


library_writer_MS2<-function(raw_data_files,metadafile,output_library,rt_search = 12,ppm_search = 20,search_adduct = T,MS2_type = c("DDA","Targeted"),baseline = 1000, normalized=T, old_library=NULL){
  
  # This function batch processes raw data files, and write them into a MGF file. It also takes into account historical library
  
  spectrum_list=list()
  metadata=c()
  NN=0 
  
  # Read from previous library:
  
  if (!is.null(old_library)){
    old_dat=readMgfData(old_library, verbose = FALSE)
    spectrum_list=Mgf2Splist(old_dat) 
    metadata=fData(old_dat)
    metadata=metadata[,1:(ncol(metadata)-1)]
    NN = length(spectrum_list)
  }
  
  # Batch process:
  
  FF = length(raw_data_files)
  
  for (ff in 1:FF){
    
    dat2 = process_MS2(raw_data_files[ff],metadafile,rt_search,ppm_search,search_adduct, MS2_type, baseline)
    nls= length(dat2$sp) # Added library size
    if (length(nls)>0){
      for (n in 1:nls){spectrum_list[[NN+n]]=dat2$sp[[n]]} # Update spectrum list
      metadata=rbind(metadata,dat2$metadata) # Update meta data
      NN = NN+nls} # Update size
  }
  
  # Filter: only keep the highest TIC
  valid=filter_scans(metadata)
  metadata = metadata[valid,]
  metadata = cbind.data.frame(metadata,SCANS=1:nrow(metadata))
  spectrum_list = spectrum_list[valid]
  
  # Save results:
  write.table(metadata,paste0(output_library,".txt"),col.names = T,row.names=F,dec=".",sep="\t")
  save(spectrum_list,metadata,file=paste0(output_library,".rData"))
  writeMGF2(spectrum_list,metadata,output_library)
}


Mgf2Splist<-function(MGFdat){
  
  # From a MSnBase object to a list of spectra m/z intensity
  N=length(MGFdat)
  spectrum_list=list()
  for (i in 1:N){spectrum_list[[i]]=cbind(MGFdat[[i]]@mz,MGFdat[[i]]@intensity)}
  return(spectrum_list)
}


filter_scans<-function(metadata){
  
  # For each feature find the scan with highest TIC, return index of filtered scans
  
  valid=c()
  ID_list=unique(metadata$ID)
  for (ID in ID_list){
    selected_rows = which(metadata$ID==ID)
    tics=metadata$TIC[selected_rows]
    wm=which.max(tics)
    valid=c(valid,selected_rows[wm])
  }
  return(valid)
}
