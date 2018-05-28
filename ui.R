#ui.R is running----
print("Run ui.R")

#init the app
source("./config.R")

fluidPage(
  theme = "style.css",
  title = 'CNViewer', #webpage title
  h3('Shiny Copy Number Variation'), #top-left title
  conditionalPanel("output.fileUploaded != true",
                   selectInput("getHgVersion", "hg build:", c('hg18', "hg19", 'hg38'), selected = "hg19", width="100px")
  ),
  fluidRow(
    column(2,
           fileInput("cnvFile", "Choose CNV File", accept = c('text/csv') ), # browse cnv file
           verbatimTextOutput("inputFileText")

    ) #cnv file input status

  ),
  hr(),
  conditionalPanel("output.fileUploaded == true",
    sidebarLayout(
      #sidebarPanel----
      sidebarPanel(width =4,
                   DT::dataTableOutput('cnvTbl'),
                   actionButton("getSNPdata", "Read in SNP data"),
                   actionButton("prevCNV", "Prev CNV"),
                   actionButton("nextCNV", "Next CNV"),
                   verbatimTextOutput("SNPdataStatus", placeholder = TRUE),
                   actionButton("addCNV", "Add CNV"),
                   actionButton("setGermCNV", "Germ CNV"),
                   actionButton("setCoveredCNV", "Covered CNV"),
                   actionButton("setFalseCNV", "False CNV"),
                   actionButton("setTrueCNV", "True CNV"),
                   splitLayout(
                     cellWidths = c(85, 80),
                     actionButton("setChr", "Set Chr"),
                     textInput("chr", "", placeholder = ""),
                     actionButton("setCN", "Set CN"),
                     textInput("copyNum", "", placeholder = ""),
                     actionButton("setNormRate", "Set NR"),
                     textInput("normRate", "", placeholder = "")
                   ),
                   #delete, set dup and clear SNP data----
                   actionButton("delCNV", "Delete CNV"),
                   #actionButton("setDupCNV", "Dup CNV"),
                   actionButton("clearSNPdata", "Clear SNPdata")
      ),
      conditionalPanel("output.SNPloaded == true",
        #mainPanel----
        mainPanel(width=8,
                  fluidRow(
                    #BAF, LRR and gene list plot----
                    column(12, plotOutput('caseBAF', height = bafLrrHeight, click = "caseBAF_click",
                                          brush = brushOpts("caseBAF_brush", delay = 500, delayType ="debounce", resetOnNew = T))),
                    column(12, plotOutput('ctrlBAF', height = bafLrrHeight, click = "ctrlBAF_click",
                                          brush = brushOpts("ctrlBAF_brush", delay = 500, delayType ="debounce", resetOnNew = T))),
                    column(12, plotOutput('caseLRR', height = bafLrrHeight, click = "caseLRR_click",
                                          brush = brushOpts("caseLRR_brush", delay = 500, delayType ="debounce", resetOnNew = T))),
                    column(12, plotOutput('ctrlLRR', height = bafLrrHeight, click = "ctrlLRR_click",
                                          brush = brushOpts("ctrlLRR_brush", delay = 500, delayType ="debounce", resetOnNew = T))),
                    column(12, plotOutput('genePlot', height = genePlotHeight , click = "genePlot_click")),
                    
                    #set position by selecting near point----
                    column(2, selectInput("getChrPos", "Chr:", c("NA", 'start', 'end'))),
                    column(2, textInput("position", "Pos:", placeholder = "", width = "200")),
                    column(1, actionButton("setStart", "Set Start")),
                    column(1, actionButton("setEnd", "Set End")),
                    #select padding folds for CNV
                    column(1, selectInput("getCnvPadX", "Padding X", c(0, 2, 4, 5, 8, 10, 20), selected = 4 )),
                    #max point number for showing BAF and LRR
                    column(1, selectInput("maxSNPnum", "maxSNP", c(1000, 2000, 5000, 8000, 10000, 20000, 50000), selected = 10000 )),
                    #zoon out----
                    column(1, actionButton("zoomOut2X", "Zoom out 2X")),
                    column(1, actionButton("zoomOut5X", "Zoom out 5X")),
                    column(1, actionButton("zoomOut10X", "Zoom out 10X"))
                  ),
                  column(2, tableOutput('geneTbl') )
                  
        )
      )
    ),
    conditionalPanel("output.SNPloaded == true",
      hr(),
      splitLayout(
            cellWidths = c(100, 200, 300, 120, 120, 120, 120, 80, 80),
            actionButton("setChrGene", "Chr/Gene"),
            textInput("chrGeneName", "", placeholder = ""),
            verbatimTextOutput("spectPos", placeholder = TRUE),
            actionButton("zoomOut2XSpect", "Zoom out 2X"),
            actionButton("zoomOut5XSpect", "Zoom out 5X"),
            actionButton("zoomOut10XSpect", "Zoom out 10X"),
            actionButton("showCNV", "show CNV"),
            actionButton("noCNV", "no CNV")
      ),
      plotOutput('spectrumPlot',
                 brush = brushOpts("spect_brush", delay = 500, delayType ="debounce", resetOnNew = T))
    )
  ),
  #allow read in cnv file again----
  tags$script('$( "#cnvFile" ).on( "click", function() { this.value = null; });')
)