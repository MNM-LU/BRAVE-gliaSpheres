sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S1_CustomArraySequences.R')" 2>&1 | tee logs/S1_CustomArraySequences_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S2_AAV-lib_libraryExtraction.R')" 2>&1 | tee logs/S2_AAV-lib_libraryExtraction_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S3_libraryIdentification.R')" 2>&1 | tee logs/S3_libraryIdentification_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S4_Fragment_translation.R')" 2>&1 | tee logs/S4_Fragment_translation_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S5_AAV-lib_tissueExtraction.R')" 2>&1 | tee logs/S5_AAV-lib_tissueExtraction_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S6_generateCompleteLibraryRanges.R')" 2>&1 | tee logs/S6_generateCompleteLibraryRanges_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S7_NormalizedGeneIdentification.R')" 2>&1 | tee logs/S7_NormalizedGeneIdentification_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S8_PlotAllGenesCoverage.R')" 2>&1 | tee logs/S8_PlotAllGenesCoverage_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S9_FullGeneHeatmap.R')" 2>&1 | tee logs/S9_FullGeneHeatmap_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S10_generateLibAnalysisPlots.R')" 2>&1 | tee logs/S10_generateLibAnalysisPlots_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S11_slidingMeanTopHits.R')" 2>&1 | tee logs/S11_slidingMeanTopHits_runLog.txt &
wait; sudo -u rstudio Rscript -e "setwd('/home/rstudio')" -e "rmarkdown::render('S12_PlotThreeSampleCoverage.R')" 2>&1 | tee logs/S12_PlotThreeSampleCoverage_runLog.txt &
