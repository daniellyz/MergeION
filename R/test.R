
#MSnbase
#tools
#gWidgetstcltk

# Create library:

rt_search = 12 # retention time tolerance
ppm_search = 10  # mass tolerance (ppm)

raw_data_files = c("NA_170405_MAS006_10.mzML","TESTMIX2_180504_MAS011_06.mzML","JNJ42165279_171214_MAS006_14.mzML")
MS2_type = c("DDA","Targeted","Targeted")
metadata_file = "library_metadata.csv"
input_library = NULL
output_library = "library_V1.mgf"
library1=library_generator(raw_data_files,metadata_file,mslevel=c(1,2),MS2_type=MS2_type,rt_search,ppm_search,search_adduct = F,baseline = 1000,normalized = T, input_library, output_library)

# Update library:

raw_data_files="GMP_R601592_150925_MAS006_04.mzXML"
MS2_type = "DDA"
metadata_file = "library_metadata_GMP.csv"
input_library = "library_V1.mgf"
output_library = "library_V2.mgf"

library2=library_generator(raw_data_files,metadata_file,mslevel=c(1,2),MS2_type=MS2_type,rt_search,ppm_search,search_adduct = F,baseline = 1000,normalized = T, input_library, output_library)
