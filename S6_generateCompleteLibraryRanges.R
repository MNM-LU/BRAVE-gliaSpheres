
#' ---
#' title: "Generate a complete library range object"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This short script generates a lookup scoring table of the AAV plasmid library so that it follows the same structure as the mRNA samples so that they can be compared for coverages.  


suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


#' Generate library rage object
#' ============================
#+ Generating library rages.......
load("data/alignedLibraries.rda")
load("data/LUTdna.rda")
load("data/multipleContfragmentsComplete.rda")
setkey(output.Table,LUTnr)
setkey(LUT.dna,LUTnr)
output.Table <- output.Table[LUT.dna,nomatch=0]
output.Table[,c("Names","i.Structure"):=NULL]
setnames(output.Table,"Sequence","fragment")
setkey(output.Table,fragment)

range.idx <- data.table(fragment=mcols(allFragments.ranges)$Sequence, 
                        idxFrag=1:length(allFragments.ranges), key="fragment")
output.Table <- output.Table[range.idx, nomatch=0, allow.cartesian=TRUE]

foundFragments.ranges <- allFragments.ranges[output.Table$idxFrag]
output.Table[,c("Reads","fragment","idxFrag","Structure","LUTnr"):=NULL]
output.Table[,RNAcount:=tCount]

mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges), output.Table)
  
saveRDS(foundFragments.ranges, file="output/completeLibraryRanges.rds")

devtools::session_info()