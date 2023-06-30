
#' ---
#' title: "Barcoded extraction and reduction from RNA samples"
#' author: "Tomas Bjorklund"
#' output:
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow identifies correct amplicons from in vivo & in vitro samples and extracts the barcode. Barcodes are then reduced using the starcode algorithm.
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(beanplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales)) #Gives the log2 ability to ggplot2
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


#' Analyze tissue RNA
#' ============================
#+ Analyzing RNA.......
strt <- Sys.time()
load("data/multipleContfragmentsComplete.rda")
load("data/alignedLibraries.rda")
load("data/LUTdna.rda")

load.list <- read.table("input/loadlist.txt", header = FALSE, skip = 0, sep="\t",
                        stringsAsFactors = FALSE, fill=TRUE)

dataDir <- "seqFiles"
colnames(load.list) <- c("Name", "BaseName","GroupName")

log.table <- data.table(Name="Name",
                        Reads=NA,
                        Purity=NA,
                        BCs=NA,
                        SCdroppedBC=NA,
                        allBCs=NA,
                        scBCs=NA)

analyzeTissue <- function(indexNr) {
#indexNr <- 1
  
  name <- unlist(strsplit(load.list$BaseName[indexNr],"/"))
  name <- name[!is.na(name)]
  if (length(name)==2){
    in.files <- list.files(paste(gsub("([\\])", "", dataDir),name[1],sep="/"), 
                           pattern=paste(name[2],"*", sep=""), full.names=TRUE)
  } else {
    in.files <- list.files(gsub("([\\])", "", dataDir), 
                           pattern=paste(name[1],"*", sep=""), full.names=TRUE)
  }
  in.files.P5 <-in.files[grep("R1",in.files)]
  in.files.P7 <- in.files[grep("R2",in.files)]

in.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
in.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
system(paste("cat '", paste(as.character(in.files.P5), collapse="' '"), "' > ", 
             in.name.P5, " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE)

system(paste("cat '", paste(as.character(in.files.P7), collapse="' '"), "' > ", 
             in.name.P7, " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE)

log.table$Name <- load.list$Name[indexNr]
name.out <- log.table$Name

# Selection of real amplicons
# ============================

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("overwrite=true k=10 rcomp=f skipr1=t qhdist=0 maskmiddle=t ",
                      "hammingdistance=1 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "CGCCACAACATCGAGGACGGCAGCGTG", sep = "") 

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) 
log.table$Purity <- strsplit(sys.out[grep("Contaminants",sys.out)],split = "\t")[[1]][2]

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7

log.table$Reads <- as.integer(system(paste("gunzip -c ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE)) #Stores the read count utilized


# Extraction of barcodes
# ============================

out.name.BC <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=12 mink=12 hammingdistance=2 ",
                        "findbestmatch=t trd=t rcomp=f skipr2=t findbestmatch=f qhdist=0 ",
                        "minavgquality=0 ordered=t maxns=0 minlength=18 maxlength=22 threads=", 
                        detectCores()," in=", shQuote(in.name.P5), " out=", out.name.BC,
                        " lliteral=", "GGCCTAGCGGCCGCTTTACTT", " rliteral=", "ATAACTTCGTATA",
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) 

log.table$BCs <- strsplit(sys.out[grep("Result:",sys.out)],split = "\t")[[1]][2]

reads.BC <- readFastq(out.name.BC)
barcodeTable <- data.table(ID=as.character(ShortRead::id(reads.BC)), 
                           BC=as.character(sread(reads.BC)), key="BC")

# Starcode based barcode reduction
# ============================

out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")

system(paste("gunzip -c ",out.name.BC," | starcode -t ",detectCores()," --print-clusters -d",
             1," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
       intern = TRUE, ignore.stdout = FALSE)

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                                     stringsAsFactors = FALSE, fill=FALSE),keep.rownames=TRUE, key="rn") 
table.BC.sc[,V2 := NULL]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3),",",fixed=TRUE), by=rn]

log.table$SCdroppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))


setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))

# Replacing barcodes with Starcode reduced versions
# ============================


setkey(table.BC.sc,BC)

barcodeTable <- barcodeTable[table.BC.sc,nomatch=0]

setnames(barcodeTable,c("BC","scBC"),c("oldBC","BC"))

setkey(barcodeTable,BC)

log.table$allBCs <- length(unique(barcodeTable$oldBC))
log.table$scBCs <- length(unique(barcodeTable$BC))

invisible(barcodeTable[,oldBC:=NULL])
setkey(output.Table,"BC")

BCcount <- data.table(as.data.frame(rev(sort(table(barcodeTable$BC))), row.names = "Var1"), keep.rownames = TRUE)
#In R versions below 3.3 remove, row.names = "Var1" to make this compatible
setnames(BCcount,colnames(BCcount),c("BC","RNAcount"))
setkey(BCcount,"BC")
foundFrags <- output.Table[BCcount,nomatch=0]
setkey(foundFrags,"LUTnr")
setkey(LUT.dna,"LUTnr")
foundFrags <- foundFrags[LUT.dna,nomatch=0]
setnames(foundFrags,"Sequence","fragment")
foundFrags[,c("Names","i.Structure"):=NULL]

matchRange <- function(idxFrag) {
  #idxFrag <- 23
  matchRanges <- which(mcols(allFragments.ranges)$Sequence == foundFrags$fragment[idxFrag])
  return(cbind(matchRanges,idxFrag))
}
match.ranges.list <- mclapply(1:nrow(foundFrags), matchRange, mc.preschedule = TRUE, 
                              mc.cores = detectCores()/2)
match.ranges <- do.call(rbind, match.ranges.list)
foundFragments.ranges <- allFragments.ranges[match.ranges[,1]]
if (ncol(match.ranges) >= 2) {
foundFrags <- foundFrags[match.ranges[,"idxFrag"],]
foundFrags[,c("Reads","fragment","Structure","LUTnr"):=NULL]
mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges),foundFrags)
o = order(-mcols(foundFragments.ranges)$RNAcount)
foundFragments.ranges <- foundFragments.ranges[o]
saveRDS(foundFragments.ranges, file=paste("output/","found.",name.out,".rds", sep=""), 
        compress = TRUE)
}
return(log.table)
}

#'Analysis summary
#'============================

all.logs <- lapply(1:nrow(load.list), analyzeTissue)
all.logs <- rbindlist(all.logs, use.names=FALSE )
knitr::kable(all.logs, format = "latex", longtable = T, booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) #Cleanup of temp files

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()