
#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.6in
#' fontsize: 9.5pt
#' ---

#' This script presents overview ploits and top candidates.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 8, fig.height = 10.5, fig.align = 'center') #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
strt1<-Sys.time()

#'Generation of infective library
#'===================
all.samples <- readRDS("data/allSamplesDataTable.RDS")

all.samples$Group[all.samples$Group== "mRNA_3cpc_HEK293T"] <- "mRNA_3cpc_HEK293T"
all.samples$Group[all.samples$Group== "mRNA_30cpc_HEK293T"] <- "mRNA_30cpc_HEK293T"

#'Plotting function
#'===================

plotPair <- function(topSample,bottomSample,size.bin=1,winWidth=1,NormalizePlot=TRUE, PlotBC=TRUE) {
  # Select samples
  #===================
  
  # topSample <- "mRNA_30cpc_Organoid_MD114_R"
  # bottomSample <- "mRNA_3000cpc_Organoid_MD101_R"
  # filterBC <- FALSE
  # filterAnimal <- FALSE
  # AnimaladjustPlot <- FALSE
  # NormalizePlot <- TRUE
  # size.bin <- 1
  # winWidth=1
  # PlotBC=TRUE
  
  fill.values <- eval(parse(text=paste("c(", topSample,"= rgb(38,64,135, maxColorValue = 255), ",
                                       bottomSample,"= rgb(157,190,217, maxColorValue = 255))",sep="")))
  setkey(all.samples,Group)
  select.samples <- all.samples[J(names(fill.values))] #Select the two compared groups
  select.samples[,RNAcount:=log2(RNAcount+1)]
  
  # if (PlotBC){
  # select.samples[,c("meanCount","SDcount","minCount"):=list(mean(RNAcount),
  #                                                sd(RNAcount),
  #                                                min(RNAcount)),by="Group"]
  # select.samples <- select.samples[RNAcount>=minCount+(2*SDcount),] #-SDcount
  # select.samples[,c("meanCount","SDcount","minCount"):=NULL]
  # }
  
  
  
  
  if (winWidth > 0) {
    setorder(select.samples,Group,GeneName,start,width)
    
    windowTable <- select.samples[,c("GeneName","start","width"), with = FALSE]
    windowTable <- unique(windowTable, by=c("GeneName","start","width"))
    windowTable <- windowTable[,(seq(width-winWidth+1)+start-1),by=c("GeneName","start","width")]
    setnames(windowTable,"V1","winStart")
    windowTable[,winEnd:=winStart+winWidth-1]
    setkeyv(windowTable,c("GeneName","start","width"))
    setkeyv(select.samples,c("GeneName","start","width"))
    select.samples.windowBin <- select.samples[windowTable, allow.cartesian=TRUE]
    select.samples.windowBin[,AAproc:=winStart/seqlength*100]
    
    setkey(select.samples.windowBin,Group)
    select.samples.windowBin <- select.samples.windowBin[J(names(fill.values))] #Select the two compared groups
    setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))
    select.samples.windowBin <- select.samples.windowBin[, list(Overlaps=.N,
                                                                seqlength=min(seqlength),
                                                                AAproc = min(AAproc),
                                                                BC = paste(t(BC), collapse=","),
                                                                Animals = paste(t(Animals), collapse=","),
                                                                LUTnrs = paste(t(LUTnrs), collapse=","),
                                                                RNAcount = sum(RNAcount)
    ), by=c("Group","GeneName","winStart","winEnd")]
    
    plot.data.dt <- unique(select.samples.windowBin, by=c("Group","GeneName","winStart","winEnd"))
    
  } else {
    plot.data.dt <- data.table::copy(select.samples)
  }
  
  #===================
  #Binning of data
  #===================
  FullLength <- 100
  position <- seq(0,FullLength,size.bin)
  plot.data.dt[,bin:=findInterval(AAproc, position)]
  
  plot.data.bin <- plot.data.dt[, list(.N,seqlength=min(seqlength),
                                       AAproc = position[findInterval(mean(AAproc),position)],
                                       BCsum=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                       AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                       LUTnrs = paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                       NormCount = sum(RNAcount)/seqlength*FullLength
  ), by=c("Group","GeneName","bin")]
  plot.data.bin <- unique(plot.data.bin, by=c("Group","GeneName","bin"))
  
  plot.data.bin[,BCanim:=as.double(BCsum*AnimalCount)]
  
  
  
  #===================
  #Filtration parameters
  #===================
  
  if (NormalizePlot) {
    plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount)
    plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount)
  }
  
  
  if (PlotBC && NormalizePlot) {
    plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim)
    plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim)
  }
  
  
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount*-1 #This line flips the values for the second group
  
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim*-1 #This line flips the values for the second group
  
  #===================
  #Output plot
  #===================
  
  if (PlotBC) {outVar <- "BCanim"} else {outVar <- "NormCount"}
  
  plot.out <- eval(parse(text=paste("ggplot(plot.data.bin,aes(x=AAproc,y=",outVar,", fill = Group, vjust=-20))",sep="")))
  plot.out <- plot.out + 
    geom_bar(stat="identity", position="identity")+theme_bw() +
    scale_fill_manual(name = "Library", values = fill.values) +
    scale_colour_manual(name = "Library", values = fill.values) +
    scale_x_continuous(limit=c(0,100), breaks=c(seq(0,100,20)), expand =c(0,0)) +
    facet_wrap(~ GeneName, ncol=5)+   
    theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),
          legend.position="bottom",
          legend.spacing=unit(0,"cm"),
          legend.key.height=unit(0, "cm"),
          plot.background=element_rect(fill="white"),
          axis.text = element_text(size = rel(0.45)),
          axis.ticks = element_line(size = rel(0.5)),
          axis.ticks.length = unit(.05, "cm"),
          strip.text.x = element_text(size = rel(0.5), colour = "black", 
                                      angle = 0, lineheight=3, vjust=-20),
          strip.background = element_blank(),
          panel.spacing.y = unit(0, "cm"),
          panel.spacing.x = unit(0, "cm"))
  
  
  #===================
  # Sort and select top samples
  #===================
  
  select.samples.binPos <<- select.samples
  setkeyv(select.samples.binPos,c("Group","structure","Sequence"))
  setorder(select.samples.binPos,Group,structure,Sequence,GeneName)
  select.samples.binPos <- unique(select.samples.binPos, by=c("Group","structure","Sequence")) 
  #Due to key, this removes replicates if identical sequence mapped to multiple genes
  
  setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))
  select.samples.binPos[,c("BCcount","NormCount","AnimalCount","LUTnrs","mainStruct","mismatches"):=
                          list(length(table(strsplit(paste(t(BC), collapse=","), ","))),
                               sum(NormCount),
                               length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                               paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                               paste(unique(structure), collapse=","),
                               median(mismatches)), by=key(select.samples.binPos)]
  
  select.samples.binPos <- unique(select.samples.binPos, by=c("Group","NormCount","LUTnrs"))
  select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount",
                                                    "BCcount","AnimalCount","LUTnrs","mainStruct",
                                                    "mismatches"), with = FALSE]
  
  if (PlotBC){
    select.samples.binPos[,BCanim:=BCcount*AnimalCount]
    setorder(select.samples.binPos,Group,-BCanim,-BCcount,-AnimalCount,-NormCount)
  } else {
    setorder(select.samples.binPos,Group,-NormCount,-BCcount,-AnimalCount)
  }
  setkey(select.samples.binPos,Group)
  select.samples.top <- select.samples.binPos[, head(.SD, 50), by=Group]
  top.sample <- select.samples.top[J(names(fill.values)[1])]
  bottom.sample <- select.samples.top[J(names(fill.values)[2])]
  top.sample[,c("Group"):=NULL]
  setnames(top.sample, "GeneName", names(fill.values)[1])
  bottom.sample[,c("Group"):=NULL]
  setnames(bottom.sample, "GeneName", names(fill.values)[2])
  
  out.list <- list(plot=plot.out,
                   plotBin=plot.data.bin,
                   top=top.sample,
                   bottom=bottom.sample)
  
  return(out.list)
}

#'Analyze samples
#'===================

#===================
# Sample plotting
#===================

#' Binning analysis version 1
#+ echo=FALSE

plotPair("mRNA_All","DNA_pscAAVlib",PlotBC=FALSE)$plot

plotPair("DNA_pscAAVlib_Prep2","DNA_pscAAVlib",PlotBC=FALSE)$plot
plotPair("DNA_AAVlib_DNAse_3cpc","DNA_AAVlib_DNAse_30cpc",PlotBC=FALSE)$plot
plotPair("DNA_AAVlib_DNAse_30cpc","DNA_pscAAVlib_Prep2",PlotBC=FALSE)$plot


out.plot.list <- plotPair("mRNA_30cpc_Trsp","mRNA_30cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
out.plot.list <- plotPair("mRNA_30cpc_Th","mRNA_30cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_30cpc_Ctx","mRNA_30cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_30cpc_SN","mRNA_30cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

plotPair("mRNA_30cpc_SN","mRNA_30cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_Ctx","mRNA_30cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_SN","mRNA_30cpc_Ctx",PlotBC=FALSE)$plot

out.plot.list <- plotPair("mRNA_3cpc_Trsp","mRNA_3cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_3cpc_Th","mRNA_3cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_3cpc_Ctx","mRNA_3cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_SN","mRNA_3cpc_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


plotPair("mRNA_3cpc_SN","mRNA_3cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Ctx","mRNA_3cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_SN","mRNA_3cpc_Ctx",PlotBC=FALSE)$plot

# 3cpc vs 30cpc analysis
plotPair("mRNA_3cpc_Trsp","mRNA_30cpc_Trsp",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Str","mRNA_30cpc_Str",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Th","mRNA_30cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Ctx","mRNA_30cpc_Ctx",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_SN","mRNA_30cpc_SN",PlotBC=FALSE)$plot

# 8wks vs 4wks analysis
plotPair("mRNA_30cpc_Str_4wks","mRNA_30cpc_Str",PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_Th_4wks","mRNA_30cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_Ctx_4wks","mRNA_30cpc_Ctx",PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_SN_4wks","mRNA_30cpc_SN",PlotBC=FALSE)$plot

plotPair("mRNA_3cpc_Str_4wks","mRNA_3cpc_Str",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Th_4wks","mRNA_3cpc_Th",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Ctx_4wks","mRNA_3cpc_Ctx",PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_SN_4wks","mRNA_3cpc_SN",PlotBC=FALSE)$plot

# In vitro analysis
out.plot.list <- plotPair("mRNA_3cpc_pNeuron","mRNA_30cpc_pNeuron",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_HEK293T","mRNA_30cpc_HEK293T",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

plotPair("mRNA_All","DNA_pscAAVlib",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Str","mRNA_30cpc_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Th","mRNA_30cpc_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Ctx","mRNA_30cpc_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_SN","mRNA_30cpc_SN",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("mRNA_30cpc_SN","mRNA_30cpc_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_Ctx","mRNA_30cpc_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_30cpc_Th","mRNA_30cpc_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("mRNA_3cpc_SN","mRNA_3cpc_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Ctx","mRNA_3cpc_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Th","mRNA_3cpc_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("mRNA_30cpc_Trsp","mRNA_30cpc_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Trsp","mRNA_3cpc_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_Trsp","mRNA_30cpc_Trsp",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("mRNA_3cpc_pNeuron","mRNA_30cpc_pNeuron",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("mRNA_3cpc_HEK293T","mRNA_30cpc_HEK293T",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

# Binning analysis version 2
plotPair("mRNA_All","DNA_pscAAVlib")$plot
plotPair("DNA_pscAAVlib_Prep2","DNA_pscAAVlib")$plot
plotPair("DNA_AAVlib_DNAse_3cpc","DNA_AAVlib_DNAse_30cpc")$plot
plotPair("DNA_AAVlib_DNAse_30cpc","DNA_pscAAVlib_Prep2")$plot

# 30cpc analysis
out.plot.list <- plotPair("mRNA_30cpc_Trsp","mRNA_30cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_30cpc_Th","mRNA_30cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_30cpc_Ctx","mRNA_30cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_30cpc_SN","mRNA_30cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

plotPair("mRNA_30cpc_SN","mRNA_30cpc_Th")$plot
plotPair("mRNA_30cpc_Ctx","mRNA_30cpc_Th")$plot
plotPair("mRNA_30cpc_SN","mRNA_30cpc_Ctx")$plot

# 3cpc analysis
out.plot.list <- plotPair("mRNA_3cpc_Trsp","mRNA_3cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_Th","mRNA_3cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_Ctx","mRNA_3cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_SN","mRNA_3cpc_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

plotPair("mRNA_3cpc_SN","mRNA_3cpc_Th")$plot
plotPair("mRNA_3cpc_Ctx","mRNA_3cpc_Th")$plot
plotPair("mRNA_3cpc_SN","mRNA_3cpc_Ctx")$plot


plotPair("mRNA_3cpc_Trsp","mRNA_30cpc_Trsp")$plot
plotPair("mRNA_3cpc_Str","mRNA_30cpc_Str")$plot
plotPair("mRNA_3cpc_Th","mRNA_30cpc_Th")$plot
plotPair("mRNA_3cpc_Ctx","mRNA_30cpc_Ctx")$plot
plotPair("mRNA_3cpc_SN","mRNA_30cpc_SN")$plot


plotPair("mRNA_30cpc_Str_4wks","mRNA_30cpc_Str")$plot
plotPair("mRNA_30cpc_Th_4wks","mRNA_30cpc_Th")$plot
plotPair("mRNA_30cpc_Ctx_4wks","mRNA_30cpc_Ctx")$plot
plotPair("mRNA_30cpc_SN_4wks","mRNA_30cpc_SN")$plot


plotPair("mRNA_3cpc_Str_4wks","mRNA_3cpc_Str")$plot
plotPair("mRNA_3cpc_Th_4wks","mRNA_3cpc_Th")$plot
plotPair("mRNA_3cpc_Ctx_4wks","mRNA_3cpc_Ctx")$plot
plotPair("mRNA_3cpc_SN_4wks","mRNA_3cpc_SN")$plot


out.plot.list <- plotPair("mRNA_3cpc_pNeuron","mRNA_30cpc_pNeuron")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

out.plot.list <- plotPair("mRNA_3cpc_HEK293T","mRNA_30cpc_HEK293T")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()
knitr::kable(out.plot.list$bottom, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()

#Novel Organoid samples

out.plot.list <- plotPair("mRNA_30cpc_Organoid_MD114","mRNA_3000cpc_Organoid_MD101")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_30cpc_Organoid_MD114","DNA_AAVlib_DNAse_30cpc")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_3000cpc_Organoid_MD101","DNA_AAVlib_DNAse_30cpc")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()



out.plot.list <- plotPair("mRNA_30cpc_Organoid_MD114","DNA_pscAAVlib")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_3000cpc_Organoid_MD101","DNA_pscAAVlib")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_30cpc_Organoid_MD114","mRNA_All")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


out.plot.list <- plotPair("mRNA_3000cpc_Organoid_MD101","mRNA_All")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>%   landscape()


#+ echo=TRUE
print("Total analysis time:")
print(Sys.time()-strt1)

devtools::session_info()

