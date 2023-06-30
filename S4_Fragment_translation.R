
#' ---
#' title: "Library fragment alignment"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow identifies correct fragments from the Cre-recombined AAV plasmid library and aligns them to the CustumArray ordered nucleotide fragments using Blastn. Consistant mutations in each fragment/barcode combination are also registered as is the putity of each barcode.


suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)

strt1<-Sys.time()

#' Load the trimmed reads
#' ============================
#+ Loading reads.......

load("data/LUTdna.rda")

fragments.file <- "data/fragments_AAVlibrary_complete.fastq.gz"
barcodes.file <- "data/barcodes_AAVlibrary_complete.fastq.gz"

reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)


#' Make CustomArray reference index for Blast
#' ============================
#+ Making Bowtie index.......

LUT.fa <- tempfile(pattern = "LUT_", tmpdir = tempdir(), fileext = ".fa")
LUT.seq = ShortRead(DNAStringSet(LUT.dna$Sequence), BStringSet(1:length(LUT.dna$LUTnr)))
writeFasta(LUT.seq,LUT.fa)

#'Save unique fragments as fasta file
#'===================
unique.reads <- unique(sread(reads.trim))

#'Select subset
#'===================
#unique.reads <- unique.reads[sample(length(unique.reads), 50000)]

unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))
fragments.unique.fa <- tempfile(pattern = "FragUnique_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(unique.reads,fragments.unique.fa)


#'Align against the library using blast
#'===================

blast.db <- tempfile(pattern = "blastDB_", tmpdir = tempdir(), fileext = ".db")
blast.out <- tempfile(pattern = "blastOut_", tmpdir = tempdir(), fileext = ".txt")

sys.out <-  system(paste("makeblastdb -in ", LUT.fa,
                         " -out ",blast.db," -dbtype nucl -title LUT -parse_seqids 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("blastn database generation")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:(nrow(sys.out)),], format =  "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")

sys.out <-  system(paste("export SHELL=/bin/sh; cat ",fragments.unique.fa," | parallel --block ",
                         floor(length(unique.reads)/detectCores()),
                         " --recstart '>' --pipe 'blastn -max_target_seqs 25 -word_size 7",
                         " -num_threads 1 -outfmt 10 -db ", blast.db,
                         " -query - '> ", blast.out, " 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) #

# table.blastn <- data.table(read.table(blast.out, header = FALSE, skip = 0, sep=";",
#                                       stringsAsFactors = FALSE, fill=FALSE) , keep.rownames=FALSE, key="V1")

system(paste("gzip -c ", blast.out, " > ./data/blastOutput.csv.gz", sep=" "))

table.blastn <- data.table(scan(file="./data/blastOutput.csv.gz", what="character", sep=";") , keep.rownames=FALSE, key="V1")


if (length(grep("Warning",table.blastn$V1)) != 0) {
  warnings.out <- unique(table.blastn[grep("Warning",table.blastn$V1),])
  table.blastn <- table.blastn[-grep("Warning",table.blastn$V1),]
  setnames(warnings.out,"V1", c("blastn Warnings"))
  # knitr::kable(warnings.out[1:(nrow(warnings.out)),], format = "markdown")
}

table.blastn[,c("Reads","Sequence","identity","alignmentLength","mismatches",
                "gapOpens", "q_start", "q_end", "s_start", "s_end", 
                "evalue","bitScore") := tstrsplit(V1,",",fixed=TRUE),]



table.blastn[,Reads:= as.character(sread(unique.reads[as.integer(Reads)]))]
table.blastn[,Sequence:= as.character(sread(LUT.seq[as.integer(Sequence)]))]
setkey(table.blastn,Sequence)
setkey(LUT.dna,Sequence)
table.blastn<- table.blastn[LUT.dna, nomatch=0]
table.blastn[,c("V1","identity","alignmentLength","gapOpens", "q_start", 
                "q_end", "s_start", "s_end", "evalue","Sequence","Names"):=NULL]
gc() #garbage collection to reduce memory foot print. Can be removed for speed

table.blastn[,bitScore:= as.numeric(bitScore)]
table.blastn[,mismatches:= as.numeric(mismatches)]

setkeyv(table.blastn,c("Reads","LUTnr"))
setorder(table.blastn,Reads,LUTnr,-bitScore) #This makes sure that a fragment is only aligned once to the reference in the top ten matches
table.blastn <- unique(table.blastn, by=c("Reads","LUTnr"))

gc() #garbage collection to reduce memory foot print. Can be removed for speed

table.blastn.topHit <- table.blastn[table.blastn[, .I[which.max(bitScore)], by="Reads"]$V1] # Select only rows with the highest bitScore

full.table <- data.table(Reads=as.character(sread(reads.trim)),
                         BC=as.character(sread(reads.BC)),
                         key="Reads")
all.reads <- nrow(full.table)

full.table <- full.table[table.blastn.topHit,nomatch=0] # Merge reads with the top hit alignment

print(paste("Alignment percentage:", percent(nrow(full.table)/all.reads)))

#' Starcode based barcode reduction
#' ============================
#+ Reducing barcodes.......

out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")

command.args <- paste("-c ",barcodes.file," | starcode -t ",detectCores()-1," --print-clusters -d",
                      1," -r5 -q -o ", out.name.BC.star, sep = "")

system2("gunzip", args=command.args, stdout=TRUE, stderr=TRUE)

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                                     stringsAsFactors = FALSE, fill=FALSE),keep.rownames=TRUE, key="rn") #, nrows = 1000
table.BC.sc[,V2 := NULL]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3),",",fixed=TRUE), by=rn]

SC.droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
print(paste("Dropped BCs in Starcode:", SC.droppedBC))

#rm(reads.BC,reads.trim)

setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))

DNA_pscAAVlib <- table.BC.sc
save(DNA_pscAAVlib, file="data/scBC_DNA_pscAAVlib.rda")


#' Replacing barcodes with Starcode reduced versions
#' ============================

setkey(full.table,BC)
setkey(table.BC.sc,BC)
full.table <- full.table[table.BC.sc,nomatch=0]
#full.table <- merge(full.table,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
#rm(table.BC.sc)

setnames(full.table,c("BC","scBC"),c("oldBC","BC"))

setkey(full.table,BC)

RetainedBC <- length(unique(full.table$oldBC))
scBC <- length(unique(full.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))


table.frag <- data.table(as.data.frame((rev(sort(table(full.table$oldBC))))[1:10], row.names = "Var1"), keep.rownames=TRUE)
#In R versions below 3.3 remove, row.names = "Var1" to make this compatible
setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
knitr::kable(table.frag, format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")

table.frag <- data.table(as.data.frame((rev(sort(table(full.table$BC))))[1:10], row.names = "Var1"), keep.rownames=TRUE)
#In R versions below 3.3 remove, row.names = "Var1" to make this compatible
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
knitr::kable(table.frag, format = "latex", booktabs = T) %>%
  kable_styling(latex_options = "striped")

invisible(full.table[,oldBC:=NULL])


#' Splitting reads into single-read and multi-read barcodes
#' ============================
#+ Splitting Reads.......


full.table <- full.table[order(full.table$BC),]
full.table[,mismatches:= as.numeric(mismatches)]

temp.table.single <- full.table[full.table[, .I[.N == 1], by="BC"]$V1]
temp.table.multi <- full.table[full.table[, .I[.N > 1], by="BC"]$V1]

temp.table.single[,c("mCount","tCount"):=1]
temp.table.single$Mode <- "Amb"
setkeyv(temp.table.multi,c("BC","LUTnr"))


temp.table.multi[,c("bitScore","mismatches" ,"tCount"):= list(mean(bitScore),
                                                              median(mismatches), .N), 
                 by=key(temp.table.multi)]

temp.table.multi$Mode <- "Def"
temp.table.multi <- unique(temp.table.multi)

print("Utilized reads.......")
print(nrow(full.table))
print("Whereof single reads.......")
print(nrow(temp.table.single))

#' Splitting multi-read barcodes into clean and chimeric
#' ============================
#+ Splitting Clean Reads.......

setkeyv(temp.table.multi,"BC")

temp.table.multi.clean <- temp.table.multi[temp.table.multi[, .I[.N == 1], by="BC"]$V1]
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[.N > 1], by="BC"]$V1]


temp.table.multi.clean[,mCount:=tCount]

print("Clean multi-read barcodes.......")
print(nrow(temp.table.multi.clean))
print("Chimeric multi-read barcodes.......")
print(length(unique(temp.table.multi$BC)))

#' Calculate consensus alignment of chimeric barcodes
#' ============================
#+ Calculation consensus reads .....

setkey(temp.table.multi,"BC")
temp.table.multi[, "mCount":=tCount]
temp.table.multi[, "tCount":=sum(tCount), by="BC"]

setkey(temp.table.multi,"Reads")
temp.table.multi[,c("LUTnr","bitScore","mismatches","Structure"):=NULL]
setkey(table.blastn,"Reads")
temp.table.multi <- temp.table.multi[table.blastn, nomatch=0, allow.cartesian=TRUE]

setkeyv(temp.table.multi,c("BC","LUTnr"))
temp.table.multi[,c("bitScore","mismatches" ,"mCount"):= list(max(bitScore),
                                                              median(mismatches), 
                                                              sum(mCount)), by=key(temp.table.multi)]
gc() #garbage collection to reduce memory foot print. Can be removed for speed
temp.table.multi <- unique(temp.table.multi, by=c("BC","LUTnr"))

setkeyv(temp.table.multi,"BC")
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[mCount == max(mCount)], 
                                                      by=key(temp.table.multi)]$V1] 
# Select only rows with the highest mCount

temp.table.multi <- temp.table.multi[temp.table.multi[, .I[which.max(bitScore)], 
                                                      by=key(temp.table.multi)]$V1] 
# Select only rows with the highest bitScore

temp.table.multi[temp.table.multi$mCount==1]$Mode <- "Amb"

print(paste("Number of barcodes with false mCount:",
            nrow(temp.table.multi[mCount > tCount])))

temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean)

print(paste("Total number of definitive Barcodes:", 
            length(grep("Def", temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", 
            length(grep("Amb", temp.table.multi.consensus$Mode))))
print(paste("Total number of single-read Barcodes:", 
            nrow(temp.table.single)))

output.Table <- rbind(temp.table.multi.consensus,temp.table.single)
save(output.Table, file="data/multipleContfragmentsComplete.rda")

print("Total analysis time:")
print(Sys.time()-strt1)
devtools::session_info()
