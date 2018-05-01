#server.R is running----
print("Run server.R")

#server main function----
function(input, output, session) {
  rv=reactiveValues()
  
  #define CNV padding folds
  rv$cnvPadding=NULL #cnvPadding * cnvLength +- CNV region
  
  #max point number for showing BAF and LRR
  rv$maxPt=NULL

  #define x range
  rv$xmin=NULL
  rv$xmax=NULL
  #read in cnv information----
  rv$cnvDF=NULL
  rv$inputFileText=NULL#indicate the status of input CNV file
  rv$sampleInfor=NULL
  observeEvent(input$cnvFile, {
    #fast read
    rawInput = fread(input$cnvFile$datapath, header = T, sep = "\t")
    if("caseID" %in% colnames(rawInput)){
      rv$inputFileText=NULL
      rv$cnvDF =  rawInput %>% 
        mutate(chr=ifelse(chr == 23, "X", ifelse(chr==24, "Y", chr)), note="False", updated="")
      rv$sampleInfor=rv$cnvDF %>% select(caseID, controlID, gender) %>% distinct(caseID, controlID, gender)
      #init SNP data
      rv$SNPdataStatus=F
      rv$SNPdata=NULL
      rv$getSNPdata=NULL
      rv$cnvItem=NULL
      gc()
    }else{
      rv$inputFileText="Incorrect input file!"
      rv$cnvDF=NULL
    }
  })
  
  output$fileUploaded <- reactive( !is.null(input$cnvFile) )
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  output$SNPloaded <- reactive( !is.null(rv$SNPdata) )
  outputOptions(output, 'SNPloaded', suspendWhenHidden=FALSE)
  
  #select columns to show in datatable plot----
  cnvDFout=reactive({
    if(!is.null(rv$cnvDF)){
      rv$cnvDF %>% select_(.dots=cnvDFoutHeader)
    }
    else{
      eval(
        parse(
          text = paste("data.frame(", paste(cnvDFoutHeader,"=NA", collapse = ", "),")" ) 
        )
      )
    }
  })
  
  #claim CNV start and end interactive to update figure immediately after update
  cnvStart=reactive({
    if(!is.null(rv$cnvItem))
      rv$cnvItem$start
  })
  cnvEnd=reactive({
    if(!is.null(rv$cnvItem))
      rv$cnvItem$end
  })
  
  #reactive values only rely on reactiveValue type, not reactive function value,
  #set xlab text----
  xlabText=reactive({
    if(!is.null(rv$cnvItem)){
      # isPaired= ifelse(sub("_.*", "", rv$cnvItem$caseID) == sub("_.*", "", rv$cnvItem$controlID), "Yes", "No")
      isPaired= rv$cnvItem$paired
      paste0("Paired: ", isPaired, ";    CN: ", rv$cnvItem$CN, ";    chr", rv$cnvItem$chr, ":", rv$cnvItem$start, "-", rv$cnvItem$end)
    }
  })
  #create cnv boundary line ggplot obj----
  cnvBoundGp=reactive({
    if(!is.null(rv$cnvItem))
      geom_vline("vline", xintercept=c(rv$cnvItem$start, rv$cnvItem$end), color="red")
  })
  #select case and control SNP data----
  rv$isUpdateSNPDF=FALSE
  observe({
    if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem)){
      rv$isUpdateSNPDF=TRUE
    }
  })
  caseSNPDF=reactive({
    if(rv$isUpdateSNPDF){
      cnvChr=rv$cnvItem$chr
      xmin=rv$xmin
      xmax=rv$xmax
      SNPDF=rv$SNPdata[[rv$cnvItem$caseID]] %>%
        .[Chr == cnvChr & Pos >= xmin & Pos <= xmax]
      #limit the max number of points
      if(nrow(SNPDF) > rv$maxPt)
        SNPDF=sample_n(SNPDF, rv$maxPt)
      return(SNPDF)
    }
  })
  ctrlSNPDF=reactive({
    if(rv$isUpdateSNPDF){
      cnvChr=rv$cnvItem$chr
      xmin=rv$xmin
      xmax=rv$xmax
      SNPDF=rv$SNPdata[[rv$cnvItem$controlID]] %>%
        .[Chr == cnvChr & Pos >= xmin & Pos <= xmax]
      #limit the max number of points
      if(nrow(SNPDF) > rv$maxPt)
        SNPDF=sample_n(SNPDF, rv$maxPt)
      return(SNPDF)
    }
  })
  #select gene list----
  cnvRegionGeneDF=reactive({
    if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem) ){
      cnvChr=rv$cnvItem$chr
      xmin=rv$xmin
      xmax=rv$xmax
      data.table(geneRefData) %>% 
        .[ chr == cnvChr & (between(start, xmin, xmax) | between(end, xmin, xmax) | (start<xmin & end>xmax) )]
    }
  })
  

  #____observe events----
  #read in SNP data----
  rv$getSNPdata=NULL
  rv$SNPdataStatus=F # no SNP data
  rv$SNPdata=NULL
  observeEvent(input$getSNPdata, {
    if(!is.null(rv$cnvDF) & !rv$SNPdataStatus){
      samples=unique(c(rv$sampleInfor$caseID, rv$sampleInfor$controlID))
      smpNum=length(samples)
      allFileExistStat=T
      withProgress({
        setProgress(message = "start reading", value = 0)
        rv$SNPdata=list()
        for(i in 1:smpNum){
          si=samples[i]
          siFile=paste0(dataDir, si, dataSuffix)
          if(!file.exists(siFile)){
            rv$getSNPdata=paste0("file ", siFile, " does not exist!")
            allFileExistStat=F
            break;
          }else{
            setProgress(message = paste0("Reading SNP data for sample: ", si), value = i/smpNum)
            #update to data.table
              rv$SNPdata[[si]]=data.table(read_tsv(siFile, col_types = "ccddd" , na = c('', 'NA', 'na', "NaN"), progress=F)) %>% set_colnames(cnvDataHeader) %>%
                .[!is.na(LRR) & !is.na(BAF)] %>%
                .[, c('LRR', 'normLRR', 'BAF', 'newPos', 'Chr'):=list(if_else(LRR < minLRR, minLRR, if_else(LRR > maxLRR, maxLRR, LRR)),
                                                                         if_else(LRR < 0, LRR/abs(minLRR), LRR),
                                                                         if_else(BAF < 0, 0, if_else(BAF > 1, 1, BAF)),
                                                                         Pos+chrPosOffset[Chr],
                                                                         factor(Chr, levels = unique(Chr))
                                                                      )]
          }
        }
      })
      if(allFileExistStat){
        rv$SNPdataStatus=T #SNP data is ready
        rv$getSNPdata="SNP data is ready!"
      }
      gc()
    }
  })
  
  
  #move to previous CNV----
  observeEvent(input$prevCNV, {
    if(!is.null(rv$cnvItem) ){
      if(input$cnvTbl_rows_selected > 1){
        cnvTblProxy = DT::dataTableProxy('cnvTbl')
        selectRows(cnvTblProxy, input$cnvTbl_rows_selected - 1)
      }
    }
  }
  )
  #move to next CNV----
  observeEvent(input$nextCNV, {
    if(!is.null(rv$cnvItem) ){
      isolate({cnvItemNum=nrow(cnvDFout()) })
      if(input$cnvTbl_rows_selected < cnvItemNum){
        cnvTblProxy = DT::dataTableProxy('cnvTbl')
        selectRows(cnvTblProxy, input$cnvTbl_rows_selected +1)
      }
      }
    }
  )
  
  
  #clear SNP data----
  observeEvent(input$clearSNPdata, {
    rv$SNPdataStatus=F
    rv$SNPdata=NULL
    rv$getSNPdata="SNP data is cleared!"
    # rv$cnvItem=NULL
    gc()
  })
  
  #get the cnv information and SNP dataframe, create ggplot obj with cnv boundary
  #observe cnv row selected----
  rv$cnvItem=NULL
  rv$cnvIdx=NULL
  observeEvent(input$cnvTbl_rows_selected, {
    rv$cnvIdx=input$cnvTbl_rows_selected
    rv$cnvItem=rv$cnvDF[rv$cnvIdx,]
    cnvLength=rv$cnvItem$end-rv$cnvItem$start+1
    rv$xmin=rv$cnvItem$start-cnvLength*rv$cnvPadding
    rv$xmax=rv$cnvItem$end+cnvLength*rv$cnvPadding
    rv$brushed=F
  })
  
  #observe case BAF/LRR brush range----
  rv$brushed=F
  observeEvent(input$caseBAF_brush, {
    rv$xmin=input$caseBAF_brush$xmin
    rv$xmax=input$caseBAF_brush$xmax
    rv$brushed=T
  })
  observeEvent(input$caseLRR_brush, {
    rv$xmin=input$caseLRR_brush$xmin
    rv$xmax=input$caseLRR_brush$xmax
    rv$brushed=T
  })
  
  #calculate xlim range----
  xLim=reactive({
    if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem)){
      xLimRange=range(caseSNPDF()$Pos)
      if(!rv$brushed){
        #add delta xlim in case cnv boundary is out of xlim
        #only when CNV row selected. not for brushed region selected
        xLimDelta=(xLimRange[2]-xLimRange[1])*0.01
        if(rv$cnvItem$start-xLimDelta < xLimRange[1]){
          xLimRange[1]=rv$cnvItem$start-xLimDelta
        }
        if(rv$cnvItem$end+xLimDelta > xLimRange[2]){
          xLimRange[2]=rv$cnvItem$end+xLimDelta
        }
      }
      return(xLimRange)
    }
  })
  
  #observe selected point----
  rv$clickedPos=NULL
  observeEvent(input$caseBAF_click, {
     if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem) ){
      nearPos=nearPoints(caseSNPDF(), input$caseBAF_click, xvar = "Pos", yvar = "BAF", maxpoints = 1, addDist = F)$Pos[1];
      if(is.na(nearPos))
        rv$clickedPos=NULL
      else
        rv$clickedPos=nearPos
      }
    })
  observeEvent(input$caseLRR_click, {
    if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem) ){
      nearPos=nearPoints(caseSNPDF(), input$caseLRR_click, xvar = "Pos", yvar = "LRR", maxpoints = 1, addDist = F)$Pos[1];
      if(is.na(nearPos))
        rv$clickedPos=NULL
      else
        rv$clickedPos=nearPos
      }
    })
  
  
  #observe clicked gene----
  geneDFout=reactive({
    if(!is.null(rv$SNPdata) & !is.null(rv$cnvItem) & !is.null(input$genePlot_click)){
      clickXpos=input$genePlot_click$x
      clickYpos=round(input$genePlot_click$y)
      cnvRegionGeneDF() %>%
        filter(start<= clickXpos, end >= clickXpos, clickYpos==lineNum) %>%
        select(-transcript, -geneColor) %>% as.data.frame() %>% select(-lineNum)
    }
  })

  #get chr start and end pos----
  observeEvent(input$getChrPos, {
      if(!is.null(rv$cnvItem) & input$getChrPos != "NA"){
        cnvChr=rv$cnvItem$chr
        chrTail=input$getChrPos
        rv$clickedPos=chrRefData %>% filter(chr== cnvChr) %>% select_(.dots=chrTail) %>% as.integer()
      }
    })
  #set CNV padding folds----
  observeEvent(input$getCnvPadX, {
      rv$cnvPadding=as.integer(input$getCnvPadX)
  })
  
  #set max SNP pt----
  observeEvent(input$maxSNPnum, {
    rv$maxPt=as.integer(input$maxSNPnum)
  })
  
  #set start and end pos to cnv----
  observeEvent(input$setStart, {
    if(!is.null(rv$cnvItem) & grepl("^\\d+$", input$position)){
       rv$cnvDF$start[rv$cnvIdx]=as.numeric(input$position)
       rv$cnvItem$start=rv$cnvDF$start[rv$cnvIdx]
       rv$cnvDF$note[rv$cnvIdx]="True"
       rv$cnvDF$updated[rv$cnvIdx]="Yes"
    }
  })
  observeEvent(input$setEnd, {
    if(!is.null(rv$cnvItem) & grepl("^\\d+$", input$position) ){
      rv$cnvDF$end[rv$cnvIdx]=as.numeric(input$position)
      rv$cnvItem$end=as.numeric(input$position)
      rv$cnvDF$note[rv$cnvIdx]="True"
      rv$cnvDF$updated[rv$cnvIdx]="Yes"
    }
  })
  
  #zoom out 2/5/10X BAF/LRR plot----
  rv$zoomRate=NULL
  observeEvent(input$zoomOut2X, {
    if(!is.null(rv$xmin)){
      rv$zoomRate=(2-1)/2
    }
  })
  observeEvent(input$zoomOut5X, {
    if(!is.null(rv$xmin)){
      rv$zoomRate=(5-1)/2
    }
  })
  observeEvent(input$zoomOut10X, {
    if(!is.null(rv$xmin)){
      rv$zoomRate=(10-1)/2
    }
  })
  observe({
    if(!is.null(rv$zoomRate) & !is.null(rv$cnvItem)){
      chrI=rv$cnvItem$chr
      chrMax=karyoList$chrLabelDF %>% filter(chr == chrI) %>% .$end
      isolate({
        xLength=rv$xmax-rv$xmin+1
        rv$xmin=rv$xmin-xLength*rv$zoomRate
        rv$xmax=rv$xmax+xLength*rv$zoomRate
        rv$xmin = ifelse(rv$xmin < 1, 1, rv$xmin)
        rv$xmax = ifelse(rv$xmax > chrMax, chrMax, rv$xmax)
        rv$zoomRate=NULL
      })
    }
  })
  
  #add cnv item----
  observeEvent(input$addCNV, {
    if(!is.null(input$cnvTbl_rows_selected) & !is.null(rv$cnvItem)){
      rv$cnvDF=rv$cnvDF %>% add_row(
        caseID=rv$cnvItem$caseID,
        controlID=rv$cnvItem$controlID,
        chr=rv$cnvItem$chr,
        start=rv$cnvItem$start,
        end=rv$cnvItem$end,
        CN=rv$cnvItem$CN,
        normRate=rv$cnvItem$normRate,
        gender=rv$cnvItem$gender,
        paired=rv$cnvItem$paired,
        .after=rv$cnvIdx)
    }
  })
  
  #delete cnv item----
  observeEvent(input$delCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF=rv$cnvDF[-input$cnvTbl_rows_selected,]
  })

  #change chr of the cnv item, init with chr start and end----
  observeEvent(input$setChr, {
    if(!is.null(input$cnvTbl_rows_selected) &  input$chr %in% chrRefData$chr){
      rv$cnvDF$chr[input$cnvTbl_rows_selected]=input$chr
      rv$cnvDF$start[input$cnvTbl_rows_selected]= chrRefData$start[chrRefData$chr == input$chr]
      rv$cnvDF$end[input$cnvTbl_rows_selected]= chrRefData$end[chrRefData$chr == input$chr]
    }
  })
  
  #set CN, normRate, true, germ and covered for cnv item----
  observeEvent(input$setCN, {
    if(!is.null(input$cnvTbl_rows_selected)){
      rv$cnvDF$CN[input$cnvTbl_rows_selected]=input$copyNum
      rv$cnvDF$note[rv$cnvIdx]="True"
      rv$cnvDF$updated[rv$cnvIdx]="Yes"
    }
  })
  observeEvent(input$setNormRate, {
    if(!is.null(input$cnvTbl_rows_selected)){
      rv$cnvDF$normRate[input$cnvTbl_rows_selected]=input$normRate
      rv$cnvDF$note[rv$cnvIdx]="True"
      rv$cnvDF$updated[rv$cnvIdx]="Yes"
    }
  })
  observeEvent(input$setTrueCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF$note[input$cnvTbl_rows_selected]="True"
  })
  observeEvent(input$setGermCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF$note[input$cnvTbl_rows_selected]="Germ"
  })
  observeEvent(input$setCoveredCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF$note[input$cnvTbl_rows_selected]="Covered"
  })
  observeEvent(input$setDupCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF$note[input$cnvTbl_rows_selected]="Dup"
  })
  observeEvent(input$setFalseCNV, {
    if(!is.null(input$cnvTbl_rows_selected))
      rv$cnvDF$note[input$cnvTbl_rows_selected]="False"
  })
  observe({
    if(!is.null(rv$cnvDF))
      write_tsv(rv$cnvDF, path = "cnv.bak")
  })
  
  #_____________________________________________________________________________________________________________________________----
  #prepare per case SNP spectrum data----
  caseSNPDFprep=reactive({
    caseID=rv$sampleInfor$caseID
    caseNum=length(caseID)
    caseSNPdata=data.table()
    for(i in 1:caseNum){
      casei=caseID[i]
      caseSNPdata=rv$SNPdata[[casei]] [, ymin:=-i] [, .(ymin, Chr, Pos, newPos, normLRR)] %>%
        bind_rows(caseSNPdata)
    }
    caseSNPdata
  })
  
  caseSNPgenomeDF=reactive({
    gc()
    rect_n(caseSNPDFprep(), maxRect, "newpos")
  })
  
  caseSNPchrDF=reactive({
    outList=list()
    chrs=unique(caseSNPDFprep()$Chr) %>% as.vector()
    chrNum=length(chrs)
    withProgress({
      setProgress(message = "Processing all chrs", value = 0)
      for(i in 1:chrNum){
        c=chrs[i]
        chrSNP=caseSNPDFprep()[Chr==c]
        chrSNP=rect_n(chrSNP, maxRect, "pos")
        outList[[c]]=chrSNP
        setProgress(message = paste0("Processing chr: ", c), value = i/chrNum)
      }
    })
    outList
  })

  #zoom out 2/5/10X/chr spectrum----
  
  observeEvent(input$zoomOut2XSpect, {
    if(rv$spectType %in% c("gene", "region") ){
      zoomRate=2
      if(rv$spectStart > 1 || rv$spectEnd < rv$chrEnd){
        spectPadding= (rv$spectEnd-rv$spectStart+1 ) * (zoomRate-1)/2
        rv$spectStart=rv$spectStart-spectPadding
        rv$spectEnd=rv$spectEnd+spectPadding
        rv$spectStart = ifelse(rv$spectStart < 1, 1, rv$spectStart)
        rv$spectEnd = ifelse(rv$spectEnd > rv$chrEnd, rv$chrEnd, rv$spectEnd)
      }
      rv$spectType = 'region'
    }
  })
  observeEvent(input$zoomOut5XSpect, {
    if(rv$spectType %in% c("gene", "region")){
      zoomRate=5
      if(rv$spectStart > 1 || rv$spectEnd < rv$chrEnd){
        spectPadding= (rv$spectEnd-rv$spectStart+1 ) * (zoomRate-1)/2
        rv$spectStart=rv$spectStart-spectPadding
        rv$spectEnd=rv$spectEnd+spectPadding
        rv$spectStart = ifelse(rv$spectStart < 1, 1, rv$spectStart)
        rv$spectEnd = ifelse(rv$spectEnd > rv$chrEnd, rv$chrEnd, rv$spectEnd)
      }
      rv$spectType = 'region'
    }
  })
  observeEvent(input$zoomOut10XSpect, {
    if(rv$spectType %in% c("gene", "region")){
      zoomRate=10
      if(rv$spectStart > 1 || rv$spectEnd < rv$chrEnd){
        spectPadding= (rv$spectEnd-rv$spectStart+1 ) * (zoomRate-1)/2
        rv$spectStart=rv$spectStart-spectPadding
        rv$spectEnd=rv$spectEnd+spectPadding
        rv$spectStart = ifelse(rv$spectStart < 1, 1, rv$spectStart)
        rv$spectEnd = ifelse(rv$spectEnd > rv$chrEnd, rv$chrEnd, rv$spectEnd)
      }
      rv$spectType = 'region'
    }
  })
  
  #observe spect brush range----
  observeEvent(input$spect_brush, {
    if(rv$spectType %in% c("chr", "gene", "region") ){
      chrOffset=chrPosOffset[rv$spectChr] %>% set_names(NULL)
      rv$spectStart=as.integer(input$spect_brush$xmin) #brush start, without offfset added
      rv$spectEnd=as.integer(input$spect_brush$xmax) #brush end, without offfset added
      rv$spectType = 'region'
    }
  })
  
  #select spectrum reigon gene list----
  spectRegionGeneDF=reactive({
    if(rv$spectType == "region" ){
      spectChr=rv$spectChr
      xmin=rv$spectStart
      xmax=rv$spectEnd
      tmp=geneRefData %>%
        filter( chr == spectChr, between(start, xmin, xmax)|between(end, xmin, xmax) | (start<xmin & end>xmax) )
      return(tmp)
    }
  })
  
  caseSpectDF=reactive({
    if(rv$SNPdataStatus){
      if(rv$spectType == "genome"){
        snpRectDF=caseSNPgenomeDF()
      }else if(rv$spectType == "chr"){
        snpRectDF=caseSNPchrDF()[[rv$spectChr]]
      }else if( rv$spectType %in% c( "region", "gene") ){#region and gene are almost the same
        caseNum=length(rv$sampleInfor$caseID)
        chr=rv$spectChr
        spectStart=rv$spectStart
        spectEnd=rv$spectEnd
        snpRectDF= caseSNPDFprep()[Chr == chr & between(Pos, spectStart, spectEnd)]
        snpNum=nrow(snpRectDF)
        if(snpNum/caseNum > maxRect*wiggleRoom){
          snpRectDF=rect_n(snpRectDF, maxRect, "pos")
        }
      }
      if(nrow(snpRectDF) >= 2){#do not plot if less than 2 points
        snpRectDF=getLRRcolor(snpRectDF)
        return(snpRectDF)
      }
    }else{
      return(NULL)
    }
  })
  
  #bottom chr/gene selected----
  rv$spectChr=NULL
  rv$spectStart=NULL #no offset
  rv$spectEnd=NULL  #no offset
  rv$geneStart=NULL #no offset
  rv$geneEnd=NULL #no offset
  rv$chrEnd=NULL
  rv$spectType="genome"
  observeEvent(input$setChrGene, {
    if(input$chrGeneName %in% chrRefData$chr){
      rv$spectChr=input$chrGeneName
      rv$spectStart=0
      rv$spectEnd=karyoList$chrLabelDF$end[karyoList$chrLabelDF$chr == rv$spectChr]
      rv$chrEnd=rv$spectEnd
      rv$spectType="chr"
    }else if(input$chrGeneName %in% geneRefData$gene){
      geneDF=geneRefData[geneRefData$gene == input$chrGeneName,]
      rv$spectChr=geneDF$chr[1]
      rv$geneStart=geneDF$start[1]
      rv$geneEnd=geneDF$end[1]
      rv$chrEnd=karyoList$chrLabelDF$end[karyoList$chrLabelDF$chr == rv$spectChr]
      geneLen=rv$geneEnd - rv$geneStart +1
      rv$spectStart=rv$geneStart-geneLen*geneAdjLengthX
      rv$spectEnd=rv$geneEnd+geneLen*geneAdjLengthX
      rv$spectStart = ifelse(rv$spectStart < 1, 1, rv$spectStart)
      rv$spectEnd = ifelse(rv$spectEnd > rv$chrEnd, rv$chrEnd, rv$spectEnd)
      rv$spectType="gene"
    }else if(grepl(pattern = "^\\w+:\\d+-\\d+$", ignore.case = T, x = input$chrGeneName)){
      inputVec=strsplit(input$chrGeneName, ":|-")[[1]]
      chr=sub(pattern = "chr", replacement = "", inputVec[1], ignore.case = T)
      start=as.integer(inputVec[2])
      end=as.integer(inputVec[3])
      if(chr %in% chrRefData$chr & start < end){
        chrMaxPos=chrLabelDF$end[chrLabelDF$chr == chr]
        if(end <= chrMaxPos){
          rv$spectChr=chr
          rv$geneStart=start
          rv$geneEnd=end
          rv$chrEnd=karyoList$chrLabelDF$end[karyoList$chrLabelDF$chr == rv$spectChr]
          geneLen=rv$geneEnd - rv$geneStart +1
          rv$spectStart=rv$geneStart-geneLen*geneAdjLengthX
          rv$spectEnd=rv$geneEnd+geneLen*geneAdjLengthX
          rv$spectStart = ifelse(rv$spectStart < 1, 1, rv$spectStart)
          rv$spectEnd = ifelse(rv$spectEnd > rv$chrEnd, rv$chrEnd, rv$spectEnd)
          rv$spectType="gene"
        }
      }
    }else{
      rv$spectType="genome"
    }
  })
  
  #___render output----
  #render cnv file status---- 
  output$inputFileText = renderText({
    if(!is.null(rv$inputFileText))
      rv$inputFileText
  })
  
  #render cnv datatable, ----
  #proxy is used to include the replaceData function to smooth the table updating
  output$cnvTbl = DT::renderDataTable( isolate(cnvDFout()),
                                       rownames = F,
                                       selection=list(mode="single"),
                                       class="compact",
                                       options = list(
                                        columnDefs = list(className = 'dt-center', targets = "_all"),
                                        pageLength = 20,
                                        lengthMenu = c(10, 20, 25, 50, 100),
                                        ordering=F,
                                        searching=F,
                                        processing = FALSE)
  )
  cnvTblProxy = DT::dataTableProxy('cnvTbl')
  observe({DT::replaceData(cnvTblProxy, cnvDFout(),
                    rownames = F,
                    resetPaging = F,
                    clearSelection = "none")
  })
  
  #render SNP data reading status----
  output$SNPdataStatus=renderText(rv$getSNPdata)
  
  #render BAF/LRR/gene plot----
  output$caseBAF = renderPlot({
    if(!is.null(rv$cnvItem) & rv$SNPdataStatus ) {
      if(nrow(caseSNPDF()) >= 2){
        par(mar=bafLrrPlotMar);
        plot(0, xaxs='i', xlim=xLim(), ylim=c(0, 1), type='n', ann=F, yaxt='n', xaxt='n')
        abline(h=c(0.25, 0.5, 0.75), col='grey', lty=2);
        points(caseSNPDF()$Pos, caseSNPDF()$BAF, pch=20, cex=ptSize)
        mtext(side = 2, text = "Case BAF", line = 2, font = 2, cex=1.5 )
        abline(v=c(rv$cnvItem$start, rv$cnvItem$end), col=cnvEdgeLineCol);
        axis(2, at=seq(0, 1, length.out = 5), tick=T, cex.axis=0.75, xpd=T);
        box();
      }
    }
  })
  output$ctrlBAF = renderPlot({
    if(!is.null(rv$cnvItem) & rv$SNPdataStatus ){
      if(nrow(caseSNPDF()) >= 2){
        par(mar=bafLrrPlotMar);
        plot(0, xaxs='i', xlim=xLim(), ylim=c(0, 1), type='n', ann=F, axes = F)
        abline(h=c(0.25, 0.5, 0.75), col='grey', lty=2);
        points(ctrlSNPDF()$Pos, ctrlSNPDF()$BAF, pch=20, cex=ptSize)
        mtext(side = 2, text = 'Control BAF', line = 2, font = 2, cex=1.5 )
        abline(v=c(rv$cnvItem$start, rv$cnvItem$end), col=cnvEdgeLineCol);
        axis(2, at=seq(0, 1, length.out = 5), tick=T, cex.axis=0.75, xpd=T);
        box();
      }
    }
  })
  output$caseLRR = renderPlot({
    if(!is.null(rv$cnvItem) & rv$SNPdataStatus ){
      if(nrow(caseSNPDF()) >= 2){
        caseSMlineDF=getSmoothLine(caseSNPDF())
        ylimRange=range(caseSNPDF()$LRR)
        par(mar=bafLrrPlotMar);
        plot(0, xaxs='i', xlim=xLim(), ylim=ylimRange, type='n', ann=F, axes = F)
        points(caseSNPDF()$Pos, caseSNPDF()$LRR, pch=20, cex=ptSize)
        mtext(side = 2, text = 'Case LRR', line = 2, font = 2, cex=1.5 )
        lines(caseSMlineDF$medPos, caseSMlineDF$medLRR, col=lrrMedLineCol);
        abline(h=0, col=lrrZeroLineCol);
        abline(v=c(rv$cnvItem$start, rv$cnvItem$end), col=cnvEdgeLineCol);
        axis(2, at=pretty(ylimRange), tick=T, cex.axis=0.75, xpd=T);
        box();
      }
    }
  })
  output$ctrlLRR = renderPlot({
    if(!is.null(rv$cnvItem) & rv$SNPdataStatus ){
      if(nrow(caseSNPDF()) >= 2){
        caseSMlineDF=getSmoothLine(ctrlSNPDF())
        ylimRange=range(ctrlSNPDF()$LRR)
        par(mar=bafLrrPlotMar);
        plot(0, xaxs='i', xlim=xLim(), ylim=ylimRange, type='n', ann=F, axes = F)
        points(ctrlSNPDF()$Pos, ctrlSNPDF()$LRR, pch=20, cex=ptSize)
        mtext(side = 2, text = 'Control LRR', line = 2, font = 2, cex=1.5 )
        lines(caseSMlineDF$medPos, caseSMlineDF$medLRR, col=lrrMedLineCol);
        abline(h=0, col=lrrZeroLineCol);
        abline(v=c(rv$cnvItem$start, rv$cnvItem$end), col=cnvEdgeLineCol);
        axis(2, at=pretty(ylimRange), tick=T, cex.axis=0.75, xpd=T);
        box();
      }
    }
  })

  output$genePlot = renderPlot({
    if(!is.null(rv$cnvItem) & rv$SNPdataStatus ){
      if(nrow(caseSNPDF()) >= 2){
        par(mar=genePlotMar);
        plot(0, xaxs='i', xlim=xLim(), ylim=c(0, 4), type='n', ann=F, axes = F)
        segments(x0=cnvRegionGeneDF()$start, x1=cnvRegionGeneDF()$end,
                 y0=cnvRegionGeneDF()$lineNum, y1=cnvRegionGeneDF()$lineNum,
                 col=cnvRegionGeneDF()$geneColor, lwd=10, lend=1)
        mtext(side = 2, text = 'Gene', line = 2, font = 2, cex=1.5 )
        mtext(side = 1, text = xlabText(), line= 1, font =2, cex=1.5)
        abline(v=c(rv$cnvItem$start, rv$cnvItem$end), col='red');
        axis(2, at=0:4, tick=T, cex.axis=0.75, xpd=T);
        box();
      }
    }
  })
  
  observe({
    if(!is.null(rv$clickedPos)){
      updateTextInput(session, "position", value = rv$clickedPos)
    }
  })
  
  #render the selected gene----
  output$geneTbl = renderTable(geneDFout(),
                               align='c'
  )
  #render bottom selected reigon----
  output$spectPos=renderText({
    if(!is.null(rv$spectChr)){
      paste0('chr', rv$spectChr, ':', rv$spectStart, '-', rv$spectEnd)
    }
  })
  #render bottom plot----
  observe({
    if(rv$SNPdataStatus){
      caseID=rv$sampleInfor$caseID
      caseNum=length(caseID)
      caseLabelDF=data.frame(caseID, y=-(1:caseNum))
      karyoHeight=1
      genderDF=data.frame(gender=rv$sampleInfor$gender ) %>%
        mutate( y1=-(1:caseNum), y2=y1+1,
                genderColor=if_else(gender == 'Male', maleCol,
                                    if_else(gender == 'Female', femaleCol, unknownGenderCol)))
      
      output$spectrumPlot = renderPlot(height = caseNum*15+250, expr={
        if(rv$SNPdataStatus & !is.null(caseSpectDF())){#duplicate but necessary
          if(rv$spectType == 'genome'){ 
            spectXmin=0
            spectXmax=genomeLen
            karyoDF=karyoList$karyoDF
            chrLabelDF=karyoList$chrLabelDF
            spectXlen=spectXmax-spectXmin
            xStep=spectXlen*btmFigStepRatio
            #plot frame
            par(mar=btmPlotMar);
            plot(0, xlim=c(spectXmin, spectXmax), ylim=c(-caseNum-5, karyoHeight), type='n', ann=F, axes = F, xaxs='i')
            segments(x0=caseSpectDF()$newPos, x1=caseSpectDF()$newPos, y0=caseSpectDF()$ymin, y1=caseSpectDF()$ymin+1, col=caseSpectDF()$LRRcolor)
            #karyotype anno rect at top
            rect(xleft=karyoDF$newStart, xright=karyoDF$newEnd, ybottom=0, ytop=karyoHeight, col=alpha(karyoDF$bandColor, alpha = 0.5), border = NA)
            #chr anno at top
            text(x = chrLabelDF$mid, y = karyoHeight+0.2, labels = chrLabelDF$chr, xpd=T, adj = c(0.5, 0))
            #case ID at left side
            text(x = spectXmin-xStep, y = -(1:caseNum)+0.5, labels = caseID, xpd=T, adj = c(1, 0.5))
            #gender rect at right tail
            rect(xleft = spectXmax, xright =spectXmax+xStep, ybottom = genderDF$y1, ytop = genderDF$y2, col = genderDF$genderColor, xpd=T)
            #first vertical line and other chr dividing line, end2 is offsetted end
            segments(x0 = c(spectXmin, chrLabelDF$end2), x1 = c(spectXmin, chrLabelDF$end2), y0 = 0, y1 = -caseNum, col=caseChrSepLineCol, xpd=T)
            #case dividing line
            segments(x0 =spectXmin, x1=spectXmax, y0=-(0:caseNum), y1=-(0:caseNum), col=caseChrSepLineCol)

          }else if(rv$spectType == 'chr'){
            spectChr=rv$spectChr
            spectXmin=rv$spectStart
            spectXmax=rv$spectEnd
            karyoDF=karyoList$karyoDF %>% filter(chr == spectChr)
            spectXlen=spectXmax-spectXmin
            xStep=spectXlen*btmFigStepRatio
            #plot frame
            par(mar=btmPlotMar);
            plot(0, xlim=c(spectXmin, spectXmax), ylim=c(-caseNum-5, karyoHeight), type='n', ann=F, axes = F, xaxs='i')
            segments(x0=caseSpectDF()$Pos, x1=caseSpectDF()$Pos, y0=caseSpectDF()$ymin, y1=caseSpectDF()$ymin+1, col=caseSpectDF()$LRRcolor)
            #karyotype anno rect at top
            rect(xleft=karyoDF$start, xright=karyoDF$end, ybottom=0, ytop=karyoHeight, col=alpha(karyoDF$bandColor, alpha = 0.5), border = NA)
            #chr anno at top
            text(x = (spectXmin+spectXmax)/2, y = karyoHeight+0.2, labels = spectChr, xpd=T, adj = c(0.5, 0))
            #case ID at left side
            text(x = spectXmin-xStep, y = -(1:caseNum)+0.5, labels = caseID, xpd=T, adj = c(1, 0.5))
            #gender rect at right tail
            rect(xleft = spectXmax, xright =spectXmax+xStep, ybottom = genderDF$y1, ytop = genderDF$y2, col = genderDF$genderColor, xpd=T)
            #first vertical line and chr dividing line, end2 is offsetted end
            segments(x0 = c(spectXmin, chrLabelDF$end2), x1 = c(spectXmin, chrLabelDF$end2), y0 = 0, y1 = -caseNum, col=caseChrSepLineCol, xpd=T)
            #case dividing line
            segments(x0 =spectXmin, x1=spectXmax, y0=-(0:caseNum), y1=-(0:caseNum), col=caseChrSepLineCol)
          }else if(rv$spectType == 'gene'){
            spectChr=rv$spectChr
            geneStart=rv$geneStart
            geneEnd=rv$geneEnd
            spectXmin=rv$spectStart
            spectXmax=rv$spectEnd
            spectXlen=spectXmax-spectXmin
            xStep=spectXlen*btmFigStepRatio
            snpPerCase=nrow(caseSpectDF())/caseNum
            #plot frame
            par(mar=btmPlotMar); 
            plot(0, xlim=c(spectXmin, spectXmax), ylim=c(-caseNum-5, karyoHeight), type='n', ann=F, axes = F, xaxs='i')
            if(snpPerCase < maxRect){
              spectSNPpad=(spectXmax-spectXmin)*spectSNPpadRatio
              padWidth=if_else(spectSNPpad < spectSNPwidthMax, spectSNPpad, spectSNPwidthMax)
              tmpDF=caseSpectDF() %>% group_by(ymin) %>% arrange(Pos) %>%
                mutate(x1=Pos - padWidth, x2=Pos + padWidth, x1=if_else(x1 < lag(x2, default = 0), (lag(x2)+x1)/2, x1),
                       x2= if_else(x2 > lead(x1, default = 3e10), lead(x1), x2))
              rect(xleft = tmpDF$x1, xright =tmpDF$x2, ybottom = tmpDF$ymin, ytop = tmpDF$ymin+1, col = tmpDF$LRRcolor, border=NA)
            }else{
              segments(x0=caseSpectDF()$Pos, x1=caseSpectDF()$Pos, y0=caseSpectDF()$ymin, y1=caseSpectDF()$ymin+1, col=caseSpectDF()$LRRcolor)
            }
            #chr anno at top
            text(x = (spectXmin+spectXmax)/2, y = 0.2, labels = spectChr, xpd=T, adj = c(0.5, 0))
            #case ID at left side
            text(x = spectXmin-xStep, y = -(1:caseNum)+0.5, labels = caseID, xpd=T, adj = c(1, 0.5))
            #gender rect at right tail
            rect(xleft = spectXmax, xright =spectXmax+xStep, ybottom = genderDF$y1, ytop = genderDF$y2, col = genderDF$genderColor, xpd=T)
            #first vertical line
            segments(x0 = spectXmin, x1 = spectXmin, y0 = 0, y1 = -caseNum, col=caseChrSepLineCol, xpd=T)
            #gene region indicator h bar
            segments(x0 = geneStart, x1 = geneEnd, y0 = -caseNum-0.3, y1 = -caseNum-0.3, col=geneAnnoCol, lwd=10, lend = 1)
            #gene region indicator v lines
            segments(x0 = c(geneStart, geneEnd), x1 = c(geneStart, geneEnd), y0 = 0, y1 = -caseNum-0.2, col=geneAnnoCol)
            #case dividing line
            segments(x0 =spectXmin, x1=spectXmax, y0=-(0:caseNum), y1=-(0:caseNum), col=caseChrSepLineCol)
          }else if(rv$spectType == 'region'){
            spectChr=rv$spectChr
            spectXmin= rv$spectStart
            spectXmax= rv$spectEnd
            spectXlen=spectXmax-spectXmin
            xStep=spectXlen*btmFigStepRatio
            snpPerCase=nrow(caseSpectDF())/caseNum
            #plot frame
            par(mar=btmPlotMar);
            plot(0, xlim=c(spectXmin, spectXmax), ylim=c(-caseNum-5, karyoHeight), type='n', ann=F, axes = F, xaxs='i')
            if(snpPerCase < maxRect){
              spectSNPpad=(spectXmax-spectXmin)*spectSNPpadRatio
              padWidth=if_else(spectSNPpad < spectSNPwidthMax, spectSNPpad, spectSNPwidthMax)
              tmpDF=caseSpectDF() %>% group_by(ymin) %>% arrange(Pos) %>%
                mutate(x1=Pos - padWidth, x2=Pos + padWidth, x1=if_else(x1 < lag(x2, default = 0), (lag(x2)+x1)/2, x1),
                       x2= if_else(x2 > lead(x1, default = 3e10), lead(x1), x2))
              rect(xleft = tmpDF$x1, xright =tmpDF$x2, ybottom = tmpDF$ymin, ytop = tmpDF$ymin+1, col = tmpDF$LRRcolor, border=NA)
            }
            else{
              segments(x0=caseSpectDF()$Pos, x1=caseSpectDF()$Pos, y0=caseSpectDF()$ymin, y1=caseSpectDF()$ymin+1, col=caseSpectDF()$LRRcolor)
            }
            #chr anno at top
            text(x = (spectXmin+spectXmax)/2, y = 0.2, labels = spectChr, xpd=T, adj = c(0.5, 0))
            #case ID at left side
            text(x = spectXmin-xStep, y = -(1:caseNum)+0.5, labels = caseID, xpd=T, adj = c(1, 0.5))
            #gender rect at right tail
            rect(xleft = spectXmax, xright =spectXmax+xStep, ybottom = genderDF$y1, ytop = genderDF$y2, col = genderDF$genderColor, xpd=T)
            #first vertical line
            segments(x0 = spectXmin, x1 = spectXmin, y0 = 0, y1 = -caseNum, col=caseChrSepLineCol, xpd=T)
            #case dividing line
            segments(x0 =spectXmin, x1=spectXmax, y0=-(0:caseNum), y1=-(0:caseNum), col=caseChrSepLineCol)
            if(nrow(spectRegionGeneDF()) < 500 & nrow(spectRegionGeneDF()) >=1){
              segments(x0=spectRegionGeneDF()$start, x1=spectRegionGeneDF()$end, xdp=NULL,
                       y0=spectRegionGeneDF()$lineNum-caseNum-4, y1=spectRegionGeneDF()$lineNum-caseNum-4,
                       col=spectRegionGeneDF()$geneColor, lwd=15, lend=1)
              text(x=spectRegionGeneDF()$start-xStep*0.25, y=spectRegionGeneDF()$lineNum-caseNum-4, labels = paste0(spectRegionGeneDF()$gene, spectRegionGeneDF()$strand), adj=c(1, 0.5), col = "dodgerblue", cex=1)
              text(x = spectXmin-xStep, y = -(caseNum+2)+0.5, labels = "Gene", xpd=T, adj = c(1, 0.5))
            }
          }
        }
      })
    }
  })
}