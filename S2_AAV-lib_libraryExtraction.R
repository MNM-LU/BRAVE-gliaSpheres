
#' ---
#' title: "Extraction of Barcodes and gene fragments"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow extracts barcodes and the gene fragments synthesized with the CustomArray using bbmap2. The fragments are then suitable for alignment to reference sequences using blastn.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)

config <- read.table("input/config.txt", header = FALSE, skip = 0, sep="\t",
                     stringsAsFactors = FALSE, fill=TRUE)
colnames(config) <- c("Parameter", "Value")

#'Sequencing files
#'===================
dataDir <- config$Value[1]
in.name.P5 <- file.path(dataDir, config$Value[2])
in.name.P7 <- file.path(dataDir, config$Value[3])
name.out <- config$Value[4]
paired.alignment <- as.logical(config$Value[5])

#'Analysis parameters
#'===================
knitr::kable(config, format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")
run.subset <- as.logical(config$Value[6])
max.cores <- as.integer(config$Value[7])
subset.count <- as.integer(config$Value[8])

strt<-Sys.time()

#' Selection of real amplicons
#' ============================
#+ Selecting real amplicons.......
# This section searches the sequencing file and only select the files with valid amplicons
out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("overwrite=true k=15 rcomp=f skipr2=t qhdist=0 maskmiddle=f",
                      " hammingdistance=2 findbestmatch=f ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "GTATGTTGTTCTGGAGCGGGAGGGTGCTATTTTGCCTAGCGATAA", sep = "")          

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE)

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7


#' Extraction of a subset
#' ============================
#+ Extracting subset.......

if (run.subset){
  suppressWarnings(sampler <- FastqSampler(gsub("([\\])", "", in.name.P5), subset.count, 
                                           readerBlockSize=1e9, ordered = TRUE)) 
  set.seed(123); tmp.P5 <- yield(sampler)
  in.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
  writeFastq(tmp.P5,in.name.P5, compress=TRUE)
  rm(tmp.P5)
  suppressWarnings(sampler <- FastqSampler(gsub("([\\])", "", in.name.P7), subset.count, 
                                           readerBlockSize=1e9, ordered = TRUE)) 
  set.seed(123); tmp.P7 <- yield(sampler)
  in.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
  writeFastq(tmp.P7,in.name.P7, compress=TRUE)
  rm(tmp.P7)
}

output.Reads <- as.integer(system(paste("gunzip -c ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE)) #Stores the read count utilized
print(paste("Utilized sequences:", output.Reads))


#' Extraction of barcodes
#' ============================
#+ Extracting barcodes.......

out.name.P5 <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=18 mink=18 hammingdistance=2 findbestmatch=t ",
                        "rcomp=f findbestmatch=f qhdist=1 minavgquality=0 maxns=0 minlength=18 ",
                        "maxlength=22 threads=", detectCores()," in=", shQuote(in.name.P5), 
                        " out=", out.name.P5," lliteral=", "GGCCTAGCGGCCGCTTTACTT",
                        " rliteral=", "ATAACTTCGTATAATGTATGC",
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) 
sys.out <- as.data.frame(sys.out)

in.name.P5 <- out.name.P5


colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")
rm(sys.out)

reads.BC <- readFastq(in.name.P5)
sread(reads.BC)
(unique.BCs <- unique(sread(reads.BC)))
output.BCs <- length(unique.BCs)
print(paste("Utilized barcodes:", output.BCs))
barcodeTable <- data.table(ID=as.character(ShortRead::id(reads.BC)), BC=as.character(sread(reads.BC)))

#' Extraction of fragments
#' ============================
#+ Extracting fragments.......


out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("overwrite=true k=18 mink=18 rcomp=f qhdist=1 maskmiddle=t",
                      " hammingdistance=2 findbestmatch=t minlength=38 maxlength=78 ordered=t ",
                      "threads=", detectCores(), " in=", in.name.P7, " out=", out.name.P7,
                      " lliteral=", "AGCAACCTCCAGAGAGGCAACG",
                      " rliteral=", "CAGACAAGCAGCTACCGCAGAT", sep = "")

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) # 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 extraction of fragments")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")

in.name.P7 <- out.name.P7

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P5_singlet <- tempfile(pattern = "P5_singlet_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7_singlet <- tempfile(pattern = "P7_singlet_", tmpdir = tempdir(), fileext = ".fastq.gz")

command.args <- paste("makepairs -c 'gzip' -f ", in.name.P5," -r ", in.name.P7,
                      " -fp ", out.name.P5, " -rp ", out.name.P7, " -fs ",
                      out.name.P5_singlet, " -rs ", out.name.P7_singlet,
                      " --stats 2>&1", sep = "")
sys.out <- system2("/usr/local/bin/pairfq", args=command.args, stdout=TRUE, stderr=TRUE)
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("pairfq pair matching")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")
rm(sys.out)


system(paste("mv ", out.name.P5, " ./data/barcodes_", name.out, ".fastq.gz", sep=""))
system(paste("mv ", out.name.P7, " ./data/fragments_", name.out, ".fastq.gz", sep=""))

unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) #Cleanup of temp files

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()

