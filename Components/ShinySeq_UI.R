#### =========== SHINYSEQ USER INTERFACE ============ ####
ui = fluidPage(
  theme = shinytheme("slate"),
  navbarPage("DESeq2 Analysis",
             
             ### First tab for uploading and displaying count data, design data, specifying DESeq Parameters ###
             tabPanel("Data Upload & Analysis Parameters", 
                      
                      ## sidebar (upload request, DESeq specifications) ##
                      sidebarPanel(
                        
                        ## Data Upload ## 
                        h4("Data Upload"),
                        div(style = "margin-top: +10px"), # reduce/increase space
                        # file browsers for raw data and info data
                        fileInput("countFile", "Upload count data (.txt)", accept = ".txt"),
                        div(style = "margin-top: -20px"), 
                        fileInput("infoFile", "Upload design data (.tsv)", accept = ".tsv"),
                        div(style = "margin-top: -15px"),
                        fileInput("gffFile", "Upload General feature format (.gff)"),
                        div(style = "margin-top: -15px"),
                        # Upload button
                        actionButton("upload", "Upload!", width = '100%',class = "btn-warning"),
                        
                        ## DESeq Design ##
                        div(style = "margin-top: +45px"),
                        h4("Experimental Design"),
                        # Dropdown menu to specify design variable for DESeq Analysis
                        selectInput("variable", "Experimental Variable:", "-"),
                        # Radio buttons to specify normalization method
                        radioButtons("normMethod", 
                                     "Normalization Method:",
                                     c("Size Factor Division", "VST")),
                        #c("Size Factor Division", "VST", "Quantile Normalization")),
                        
                        # Button to start DESeq
                        actionButton("analyze", "Run DESeq!", width = '100%',class = "btn-warning"),
                        
                        ## Contrast comparison of experimental variables ## 
                        div(style = "margin-top: +45px"),
                        h4("Experimental Comparison"),
                        selectInput("contrast1", 'First variable ("baseline"):', "-"),
                        selectInput("contrast2", 'Second variable:', "-"),
                        # Slider tospecify significance level
                        sliderInput("alpha", 
                                    "Significance Level:", 
                                    min = 0.01, 
                                    max = 1, 
                                    step = 0.01, 
                                    value = 0.05),
                        actionButton("results", "Get Results!", width = '100%',class = "btn-warning"),
                        
                      ), # side bar close
                      
                      ## Main panel displaying data, results, plots ##
                      mainPanel(tabsetPanel(type = "tabs",
                                            # Raw data and info table
                                            tabPanel("Raw counts", downloadButton("downloadCounts", "Download"), tableOutput("countTable")),
                                            tabPanel("Design", tableOutput("designTable")),
                                            # Normalized data and results
                                            #tabPanel("Normalized Counts", tableOutput("normalizedTable")),
                                            tabPanel("Normalized Counts",
                                                     tabsetPanel(type = "tabs",
                                                                 tabPanel("Specified method", tableOutput("normalizedTable")),
                                                                 tabPanel("TPM-normalized", tableOutput("tpmTable"))
                                                     )
                                            ),
                                            tabPanel("Results",
                                                     tabsetPanel(type = "tabs",
                                                                 tabPanel("All", textOutput("resText_all"), downloadButton("downloadResults", "Download"), tableOutput("resTable")),
                                                                 tabPanel("Significant", textOutput("resText_sig"), downloadButton("downloadSignificant", "Download"), tableOutput("significantTable"))
                                                        )
                                                     )
                      )
                      ) # main panel close
             ), # tabPanel "Data Upload & Analysis Parameters" close
             tabPanel("Plots",
                      tabsetPanel(type = "tabs",
                                  # BOXPLOTS 
                                  tabPanel("Boxplots", plotOutput("boxplot")),
                                  # PCA
                                  tabPanel("PCA", 
                                           sidebarPanel(
                                             h4("Select up to two conditions:"),
                                             # 1 dropdown menu per condition:
                                             selectInput("pca1", "First Condition (required):", "-"),
                                             selectInput("pca2", "Second Condition (optional):", "-"),
                                             actionButton("pcaPlot", "Refresh Plot!", width = '100%', class = "btn-warning")
                                           ),
                                           mainPanel(plotOutput("pca", brush = "pca_brush"), verbatimTextOutput("pca_info"))), # interactive PCA plot
                                  # HEATMAPS
                                  tabPanel("Heatmaps", 
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Experiments", plotOutput("heatExp")), # distance heatmap of experiment data
                                                       tabPanel("Genes (Variance)",                    # tabPanel for heatmap of high variance genes
                                                                sidebarPanel(
                                                                  h4("Select amount of Genes:"),
                                                                  sliderInput("geneHeatNo",            # slider to select amount of genes
                                                                              "Number", 
                                                                              min = 5, 
                                                                              max = 200, 
                                                                              step = 1, 
                                                                              value = 20),
                                                                  h4("Plot height:"),
                                                                  sliderInput("geneHeatHeight",
                                                                              "Pixels",
                                                                              min = 400,
                                                                              max = 2000,
                                                                              step = 10,
                                                                              value = 500),
                                                                  actionButton("plotGeneHeat", "Refresh Plot!", width = '100%', class = "btn-warning")
                                                                ),
                                                                mainPanel(plotOutput("heatGene"))
                                                       ) # tabPanel "Genes" close
                                           ) # tabsetPanel close
                                  ), # tabPanel "Heatmaps" close
                                  
                                  # VOLCANO PLOT
                                  tabPanel("Volcano", 
                                           sidebarPanel(
                                             h4("Select LogFC Threshold:"),
                                             sliderInput("volcFcThr",       # slider to select threshold for the logFC of the volcanoplot
                                                         "Threshold (absolute)",
                                                         min = 0,
                                                         max = 3,
                                                         step = 0.1,
                                                         value = 1),
                                             actionButton("volcPlot", "Refresh Plot!", width = '100%', class = "btn-warning")
                                           ),
                                           mainPanel(plotOutput("volc", brush = "volc_brush"), verbatimTextOutput("volc_info")))
                      )
             ), # tabPanel "Plots" close
             ### OVERVIEW TABLE OVER UP- AND DOWNREGULATED GENES ### 
             tabPanel("Up-/Downregulation Overview",
                      sidebarPanel(
                        h4("Specify Contrast"),
                        selectInput("contrastUpDown_1", "First Experimental Group:", "-"),
                        selectInput("contrastUpDown_2", "Second Experimental Group:", "-"),
                        actionButton("addToOverview", "Add to table!", width = '100%',class = "btn-warning"),
                        actionButton("clearOverview", "Clear table!", width = '100%',class = "btn-warning"),
                      ), # sidebarPanel close
                      mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel("Table", textOutput("overviewInfo"), downloadButton("downloadOverview", "Download"), tableOutput("overviewTable")),
                                    tabPanel("Venn Diagram", plotOutput("vennDiagram")),
                                    tabPanel("UpSet Plot", plotOutput("upsetPlot"))
                          ) # tabsetPanel of mainPanel close
                        ) # mainPanel close
             ) # tabPanel "Overview" close
  ) # navBarPage close
) # ui close