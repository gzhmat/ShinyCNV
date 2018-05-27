print("Run config.R")

#install relied packages
source("misc/install.R")
#init input and library
rm(list=ls())
library(shiny)
library(magrittr)
library(scales)
library(dplyr)
library(data.table)
library(DT)
library(readr)

#change if necessary
dataDir="./SNPdata/"
dataSuffix=".txt.gz"

#global variables

#max point number for showing the spectrum/rect
maxRect=2500 #base on the screen resolution

#control the width of SNP in spectrum when expanded
spectSNPpadRatio=1/400
spectSNPwidthMax=2000

#SNP number less than wiggleRoom*maxRect will not be sampled
wiggleRoom=1.5
bafLrrHeight=140
genePlotHeight=100
# cnvPadding=3 

maxLRR=1
minLRR=-4
cnvDataHeader=c("Name", "Chr", "Pos", "LRR", "BAF")

cnvDFoutHeader=c('index', 'caseID', 'chr', 'start', 'end', 'CN', 'normRate', 'note', 'updated')

#BAF/LRR plot
ptSize=0.5
cnvEdgeLineCol='red'
lrrMedLineCol='red'
lrrZeroLineCol='skyblue'

bafLrrPlotMar=c(0.25, 4, 0.25, 0)
genePlotMar=c(3, 4, 0.25, 0)

#bottom plot----

maleCol='skyblue'
femaleCol='pink'
karyoHeight=1
unknownGenderCol='grey80'
btmPlotMar=c(2, 10, 3, 2)
caseChrSepLineCol="grey25"
geneAnnoCol="green4"
btmFigStepRatio=1/200
#define gene/region padding folds
geneAdjLengthX=10

#cytoband width for label
cytoLabelWidth=0.02
cytoLabelCex=0.75

#CNV color
CNV_loss1_col="dodgerblue"
CNV_loss2_col="blue4"
CNV_gain1_col="magenta3"
CNV_gain2_col="red4"
CNV_LOH_col="green4"
CNV_select_col="orangered"
CNV_hight=0.25

#prepare smooth line data----
getSmoothLine <- function (snpData, smoothNum=10){
  rowNum=nrow(snpData)
  snpData %>% arrange(Pos) %>% mutate(smtGroup=as.integer(1:rowNum/smoothNum)) %>%
    group_by(smtGroup) %>% summarise(medPos=median(Pos), medLRR=median(LRR) ) %>% ungroup()
}

rect_n <- function(snpData, n, type){
  pos=snpData$Pos
  if(type == "newpos"){
    pos=snpData$newPos
  }
  start=min(pos)
  end=max(pos)
  len=end-start+1
  interval=len/(n*wiggleRoom)
  #data.table
  if(type == "newpos"){
     data.table(snpData)[, label:=as.integer((newPos - start)/interval)] %>%
     .[, .(Pos=mean(Pos), newPos=mean(newPos), normLRR=median(normLRR) ), by=list(ymin, label)]
  }else{
    data.table(snpData)[, label:=as.integer((Pos -    start)/interval)] %>%
      .[, .(Pos=mean(Pos), newPos=mean(newPos), normLRR=median(normLRR)), by=list(ymin, label)]
  }
  #dplyr is too slow
  # if(type == "newpos"){
  #   snpData %>% mutate(label=as.integer((newPos - start)/interval)) %>% group_by(ymin, label) %>%
  #     summarise(Pos=mean(Pos), newPos=mean(newPos), normLRR=median(normLRR))
  # }else{
  #   snpData %>% mutate(label=as.integer((Pos - start)/interval)) %>% group_by(ymin, label) %>%
  #     summarise(Pos=mean(Pos), newPos=mean(newPos), normLRR=median(normLRR))
  # }
}

#parallel version, but slower then normal version
getCaseSNPchrDF <- function(caseSNPDFprep){
  chrList=list()
  outList=list()
  chrs=unique(caseSNPDFprep$Chr) %>% as.vector()
  chrNum=length(chrs)
  # Calculate the number of cores
  no_cores=detectCores()- 2
  myCluster = makeCluster(no_cores)
  clusterExport(myCluster, list("rect_n", "maxRect", "wiggleRoom", "%>%", "mutate", "group_by", "summarise"))
  for(i in 1:chrNum){ #do not use chrs to make the order is correct
    c=chrs[i]
    chrSNP=caseSNPDFprep %>% filter(Chr==c)
    chrList[[c]]=chrSNP
  }
  outList=parLapply(myCluster, chrList,
                    function(x){
                      rect_n(x, maxRect, "pos")
                    }
  )
  stopCluster(myCluster)
  outList
}

#read in cytoband infor----
getKaryoData <- function(file){
  karyoColors=c(acen='pink4',gneg='grey20',gpos100='grey35',gpos75='grey45',
                gpos50='grey55',gpos25='grey65',gvar='grey75',stalk='grey10')
  karyoData=read_tsv(file=file, col_types = "ccddc", progress = F) %>% mutate(chr=sub("chr", "",chr) )
  chrLen=karyoData %>% select(chr, end) %>% group_by(chr) %>% filter(end == max(end)) %>% ungroup()
  chrPosOffset= chrLen %>% mutate(posOffset=cumsum(as.double(end)), posOffset=lag(posOffset, default = 0)) %>%
    .$posOffset %>% as.vector() %>% set_names(unique(chrLen$chr))
  chrLabelDF=chrLen %>% mutate(mid=end/2+chrPosOffset[chr], end2=end+chrPosOffset[chr]) %>% as.data.frame()
  karyoDF=karyoData %>% mutate(newStart=start+chrPosOffset[chr], newEnd=end+chrPosOffset[chr], bandColor=karyoColors[stain])
  return(list(karyoDF=karyoDF, chrLabelDF=chrLabelDF, chrPosOffset=chrPosOffset))
}

#map input value to most adjacent std value
toStdVal<-function(inputVal, stdVal){
  Map(function(x) stdVal[which.min(abs(stdVal-x))], inputVal) %>% unlist()
}

#prepare LRR color data----
getLRRcolor<- function (snpData){
  upVal=round(1-log2(seq(4, 1, length.out = 101))/2, digits = 3)
  upColor=colorRampPalette(c("grey95","red4"))(101) %>% set_names(upVal)
  downVal=round(log2(seq(1, 16, length.out = 101))/4-1, digits = 3)
  downColor=colorRampPalette(c("blue4", "grey95"))(101) %>% set_names(downVal)
  snpData %>%
    mutate(normLRR=round(normLRR, digits = 3),
           LRRcolor=if_else(normLRR > 0, upColor[as.character(toStdVal(normLRR, upVal))],
                            downColor[as.character(toStdVal(normLRR, downVal))]))
}
