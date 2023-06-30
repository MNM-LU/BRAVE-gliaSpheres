# Generation of the sliding window amino acid fragments for CustomArray synthesis
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S1_CustomArraySequences.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S1_CustomArraySequences_runLog.txt")
system2("mv", args="S1_CustomArraySequences.pdf output/")

# Extraction of barcodes and fragments from the Cre-recombined plasmin library
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S2_AAV-lib_libraryExtraction.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S2_AAV-lib_libraryExtraction_runLog.txt")
system2("mv", args="S2_AAV-lib_libraryExtraction.pdf output/")

# Identification and alignments of CustomArray generated fragments to the reference sequences 
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S3_libraryIdentification.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S3_libraryIdentification_runLog.txt")
system2("mv", args="S3_libraryIdentification.pdf output/")

# S4_Fragment_translation
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S4_Fragment_translation.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S4_Fragment_translation_runLog.txt")
system2("mv", args="S4_Fragment_translation.pdf output/")

# S5_AAV-lib_tissueExtraction
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S5_AAV-lib_tissueExtraction.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S5_AAV-lib_tissueExtraction_runLog.txt")
system2("mv", args="S5_AAV-lib_tissueExtraction.pdf output/")

# S6_generateCompleteLibraryRanges
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S6_generateCompleteLibraryRanges.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S6_generateCompleteLibraryRanges_runLog.txt")
system2("mv", args="S6_generateCompleteLibraryRanges.pdf output/")

# S7_NormalizedGeneIdentification
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S7_NormalizedGeneIdentification.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S7_NormalizedGeneIdentification_runLog.txt")
system2("mv", args="S7_NormalizedGeneIdentification.pdf output/")

# S8_PlotAllGenesCoverage
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S8_PlotAllGenesCoverage.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S8_PlotAllGenesCoverage_runLog.txt")
system2("mv", args="S8_PlotAllGenesCoverage.pdf output/")

# S9_FullGeneHeatmap
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S9_FullGeneHeatmap.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S9_FullGeneHeatmap_runLog.txt")
system2("mv", args="S9_FullGeneHeatmap.pdf output/")

# S10_generateLibAnalysisPlots
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S10_generateLibAnalysisPlots.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S10_generateLibAnalysisPlots_runLog.txt")
system2("mv", args="S10_generateLibAnalysisPlots.pdf output/")

# S11_slidingMeanTopHits
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S11_slidingMeanTopHits.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S11_slidingMeanTopHits_runLog.txt")
system2("mv", args="S11_slidingMeanTopHits.pdf output/")

# S12_PlotThreeSampleCoverage
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S12_PlotThreeSampleCoverage.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S12_PlotThreeSampleCoverage_runLog.txt")
system2("mv", args="S12_PlotThreeSampleCoverage.pdf output/")

# S13_Hammock_HMM_Clustering
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S13_Hammock_HMM_Clustering.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S13_Hammock_HMM_Clustering_runLog.txt")
system2("mv", args="S13_Hammock_HMM_Clustering.pdf output/")

# S14_HMM_Selected_Clustering
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S14_HMM_Selected_Clustering.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S14_HMM_Selected_Clustering_runLog.txt")
system2("mv", args="S14_HMM_Selected_Clustering.pdf output/")

# S15_HammockHeatmap
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S15_HammockHeatmap.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S15_HammockHeatmap_runLog.txt")
system2("mv", args="S15_HammockHeatmap.pdf output/")

# S16_Hammock_HMM_Organoids
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S16_Hammock_HMM_Organoids.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S16_Hammock_HMM_Organoids_runLog.txt")
system2("mv", args="S16_Hammock_HMM_Organoids.pdf output/")

# S17_Hammock_HMM_Clustering_w-organoids
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S17_Hammock_HMM_Clustering_w-organoids.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S17_Hammock_HMM_Clustering_w-organoids_runLog.txt")
system2("mv", args="S17_Hammock_HMM_Clustering_w-organoids.pdf output/")

# S18_HMM_Selected_Clustering_w-organoids
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S18_HMM_Selected_Clustering_w-organoids.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S18_HMM_Selected_Clustering_w-organoids_runLog.txt")
system2("mv", args="S18_HMM_Selected_Clustering_w-organoids.pdf output/")

# S19_HammockHeatmap_w-organoids
sys.args <- paste("-e",shQuote("setwd('/home/rstudio')"),"-e",shQuote("rmarkdown::render('S19_HammockHeatmap_w-organoids.R')"), sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/S19_HammockHeatmap_w-organoids_runLog.txt")
system2("mv", args="S19_HammockHeatmap_w-organoids.pdf output/")

