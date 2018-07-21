# MergeION

In tandem mass spectrometry-based metabolomics, automated structure identification is usually performed by spectral library search and in silico fragmentation (smart algorithms).
The missing steps towards complete automation are the pre-processing and format conversion of raw chromatograms. 
The package fills these gaps and enables batch-processing of targeted MS/MS scans on commonly used structure elucidation software, including CSI:FingerID, MSFinder, GNPS and Metfrag. 
To achieve this, MS/MS scans from one or multiple raw chromatogram files are extracted according to m/z (and retention time) provided by users. 
They are then merged into a single "spectral library" that also contain user-provided metadata. The spectral library can be exported to different formats compatible with different software tools.
