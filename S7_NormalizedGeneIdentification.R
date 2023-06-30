
#' ---
#' title: "Normalize Library counts"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow normalizes read counts between samples to compensate for variable read depth.  

suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(multicore))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


#' Generate load list and grouping names
#' ============================
#+ Generating load list.......
strt <- Sys.time()
in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
load.list <- read.table("input/loadlist.txt", header = FALSE, skip = 0, sep="\t",
                        stringsAsFactors = FALSE, fill=TRUE)
colnames(load.list) <- c("Name", "BaseName","GroupName")
load.list <- rbind(load.list,c("completeLibraryRanges","","DNA_pscAAVlib"))
load.list <- load.list[!grepl("Untreat",load.list$Name),]

select.Cases <- c(unlist(sapply(load.list$Name, function(x) grep(x,in.names.all), simplify = TRUE)))

(in.names.all <- in.names.all[select.Cases])


grouping <- data.frame(Sample=gsub("-","_",gsub("found.|(output/)|(.rds)", "", in.names.all)),
                       Group=load.list[match(names(select.Cases),load.list$Name),"GroupName"],
                       stringsAsFactors = FALSE)

#' Load the desired alignment files and annotating group
#' ============================
#+ Loading alignments.......

loadRDS <- function(in.name) {
  #in.name <- in.names.all[42]
  this.sample <- readRDS(in.name)
  this.name <- gsub("-","_",gsub("found.|(output/)|(.rds)", "", in.name))
  this.group <- grouping[match(this.name,grouping$Sample),"Group"]
  mcols(this.sample) <- cbind(mcols(this.sample),
                              data.frame(Sample = this.name, Group=this.group,
                                         stringsAsFactors = FALSE))
  
  return(this.sample)
}

all.samples <- lapply(in.names.all, loadRDS)

all.samples <- do.call(GAlignmentsList,unlist(all.samples))
all.samples <- cbind(unlist(all.samples))[[1]]

names(all.samples) <- make.names(names(all.samples), unique=TRUE)
length.Table <- data.table(seqnames=names(seqlengths(all.samples)),
                           seqlength=seqlengths(all.samples), key="seqnames")
all.samples <- data.table(as.data.frame(all.samples), key="seqnames")
all.samples[,c("strand","qwidth","cigar","njunc","end"):=NULL]
all.samples <- all.samples[length.Table] #A data.table merge to match seqlengths to their respective seqnames
all.samples[, c("Category", "Protein", "Origin", 
                "Extra", "Number","GeneName") := tstrsplit(seqnames, ",", fixed=TRUE)]
all.samples[, c("seqnames","Protein", "Origin", 
                "Extra", "Number") := NULL]
all.samples[, GeneName := gsub("/|_","-",GeneName)]


#' Normalizing read counts to correct for variable read depth
#' ============================
#+ Normalizing RNA counts.......
setkey(all.samples,Group)
all.samples <- all.samples[RNAcount>1,] #Filters out single count reads
readCounts <- all.samples[,list(GroupCount=sum(RNAcount)), by="Group"]
readCounts[,GroupCount:=GroupCount/max(GroupCount)]
setkey(readCounts,Group)
all.samples <- all.samples[readCounts] #Merge with normalizing factor
all.samples[,RNAcount:=RNAcount/GroupCount]
setkey(all.samples,Mode)
all.samples <- all.samples["Def"]

setkey(all.samples,Group)
total.AAV.samples <- all.samples[Group!="DNA_pscAAVlib" & Group!="DNA_pscAAVlib_Prep2" & Group!="DNA_AAVlib_DNAse_3cpc" & Group!="DNA_AAVlib_DNAse_30cpc"]
#total.AAV.samples <- total.AAV.samples[!grepl("4wks",total.AAV.samples$Group)]
transported.AAV.samples.30cpc <- total.AAV.samples[grepl("mRNA_30cpc_SN|mRNA_30cpc_Th|mRNA_30cpc_Ctx",total.AAV.samples$Group)]
transported.AAV.samples.3cpc <- total.AAV.samples[grepl("mRNA_3cpc_SN|30cpc_Th|mRNA_3cpc_Ctx",total.AAV.samples$Group)]
total.AAV.samples[,Group := "mRNA_All"]
transported.AAV.samples.30cpc[,Group := "mRNA_30cpc_Trsp"]
transported.AAV.samples.3cpc[,Group := "mRNA_3cpc_Trsp"]

all.samples <- rbind(all.samples,total.AAV.samples,transported.AAV.samples.30cpc,transported.AAV.samples.3cpc)

rm(total.AAV.samples,transported.AAV.samples.30cpc,transported.AAV.samples.3cpc)

setkeyv(all.samples,c("Group","Category","GeneName","structure","start","width","Sequence","seqlength"))

all.samples <- all.samples[,j=list(bitScore=sum(bitScore*tCount)/sum(tCount),
                                   mismatches=median(mismatches),
                                   mCount=sum(mCount),
                                   tCount=sum(tCount),
                                   BC=paste(unique(BC), collapse = ","),
                                   Animals=paste(unique(Sample), collapse = ","),
                                   LUTnrs=paste(unique(LUTnr), collapse = ","),
                                   RNAcount=sum(RNAcount),
                                   NormCount=log2(sum(RNAcount)+1)*.N),
                           by=c("Group","Category","GeneName","structure","start","width","Sequence","seqlength")]

all.samples[,start:=floor((start+2)/3)]
all.samples[,width:=ceiling((width)/3)]
all.samples[,seqlength:=ceiling(seqlength/3)]
all.samples[,AA:=floor(start+(width/2))]
all.samples[,AAproc:=AA/seqlength*100]

#' Remove overhangs on the sequence based on the Structure annotation
#' ============================
#+ Removing overhangs.......

all.samples[structure == "14aa", "Sequence" := substr(Sequence,3,44)]
all.samples[structure == "22aa", "Sequence" := substr(Sequence,3,68)]
all.samples[structure == "14aaG4S", "Sequence" := substr(Sequence,15,56)]
all.samples[structure == "14aaA5", "Sequence" := substr(Sequence,15,56)]

#Change the default behavior to induce start codons and Methionine
GENETIC_CODE_ALT <- GENETIC_CODE
attr(GENETIC_CODE_ALT, "alt_init_codons") <- c("TAA","TAG")

all.samples[, Peptide := mclapply(Sequence, function(x) as.character(Biostrings::translate(DNAString(x), genetic.code=GENETIC_CODE_ALT, if.fuzzy.codon="solve")), mc.cores = detectCores())]
all.samples[,Peptide:= as.character(Peptide),]
saveRDS(all.samples, file="data/allSamplesDataTable.RDS")

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()

