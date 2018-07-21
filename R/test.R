
#MSnbase
#tools
#gWidgetstcltk

source("library_writer.r")

rt_search = 12 # retention time tolerance 
ppm_search = 10  # mass tolerance (ppm)

### Writing first file (DDA mode):

#raw_data_files = c(list.files(pattern="\\.mzML$"),list.files(pattern="\\.mzXML$"))

raw_data_files = c("NA_170405_MAS006_10.mzML","TESTMIX2_180504_MAS011_06.mzML")
metadafile = "library_metadata.csv"

# MS1:
output_library = "library_MS1_V1.mgf"
library_writer_MS1(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = T,baseline = 1000,normalized = T)

# MS2:
output_library = "library_MS2_V1.mgf"
library_writer_MS2(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = F,MS2_type = "DDA",baseline = 1000,normalized = T)

### Writing second file:

raw_data_files = c("JNJ42165279_171214_MAS006_14.mzML")
metadafile = "library_metadata.csv"

# MS1:
input_library = "library_MS1_V1.mgf"
output_library = "library_MS1_V2.mgf"
library_writer_MS1(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = F,baseline = 1000,normalized = T, old_library = input_library)

# MS2:
input_library = "library_MS2_V1.mgf"
output_library = "library_MS2_V2.mgf"
library_writer_MS2(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = F,MS2_type = "Targeted",baseline = 1000,normalized = T, old_library = input_library)

### Writing third file:

raw_data_files = c("GMP_R601592_150925_MAS006_04.mzXML")
metadafile = "library_metadata_GMP.csv"

# MS1:
input_library = "library_MS1_V2.mgf"
output_library = "library_MS1_V3.mgf"
library_writer_MS1(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = F,baseline = 1000,normalized = T, old_library = input_library)

# MS2:
input_library = "library_MS2_V2.mgf"
output_library = "library_MS2_V3.mgf"
library_writer_MS2(raw_data_files,metadafile,output_library,rt_search,ppm_search,search_adduct = F,MS2_type = "Targeted",baseline = 1000,normalized = T, old_library = input_library)
