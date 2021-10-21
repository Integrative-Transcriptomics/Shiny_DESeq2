#### =========== SHINYSEQ USER INTERFACE ============ ####
ui = fluidPage(
  tags$head(tags$style(HTML('.modal-lg {width: 100%;}'))),
  theme = shinytheme("cerulean"),
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
                        # Slider to specify significance level
                        sliderInput("alpha", 
                                    "Significance Level:", 
                                    min = 0.01, 
                                    max = 1, 
                                    step = 0.01, 
                                    value = 0.05),
                        # Button to start DESeq
                        actionButton("analyze", "Run DESeq!", width = '100%',class = "btn-warning"),
                      ), # side bar close
                      
                      ## Main panel displaying data, results, plots ##
                      mainPanel(tabsetPanel(type = "tabs",
                                            # Raw data and info table
                                            tabPanel("Raw counts", downloadButton("downloadCounts", "Download"), dataTableOutput("countTable")),
                                            tabPanel("Design", downloadButton("downloadDesign", "Download"), dataTableOutput("designTable"))
                      )
                      ) # main panel close
             ), # tabPanel "Data Upload & Analysis Parameters" close
             tabPanel("Normalization",
                      tabsetPanel(type = "tabs",
                                  # Normalized data
                                  tabPanel("Normalized Counts",
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Specified method", downloadButton("downloadNormalizedCounts", "Download"), dataTableOutput("normalizedTable")),
                                                       tabPanel("TPM-normalized", downloadButton("downloadTPMCounts", "Download"), dataTableOutput("tpmTable"))
                                           )
                                  ),
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
                                                       tabPanel("Experiments", 
                                                                sidebarPanel(h4("Size:"),
                                                                             sliderInput("experimentHeatHeight", 
                                                                                         "Height (px)",
                                                                                         min = 500,
                                                                                         max = 1500,
                                                                                         step = 10,
                                                                                         value = 750),
                                                                             sliderInput("experimentHeatWidth", 
                                                                                         "Width (px)",
                                                                                         min = 500,
                                                                                         max = 1500,
                                                                                         step = 10,
                                                                                         value = 750),
                                                                             actionButton("plotExperimentHeat", "Refresh Plot!", width = '100%', class = "btn-warning")
                                                                             ),
                                                                mainPanel(plotOutput("heatExp"))
                                                                ), 
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
                                  tabPanel("Profile Plots",
                                           sidebarPanel(
                                             h4("Add genes to profile plot"),
                                             selectizeInput("profileGene", "Gene name(s)", choices = "-", multiple = TRUE),
                                             radioButtons("averageReplicates", "Average replicates?", choices = c("Yes", "No")),
                                             radioButtons("profileErrorbars", "Add errorbars?", choices = c("Yes", "No")),
                                             actionButton("profilePlotButton", "Refresh Plot!", width = '100%', class = "btn-warning"),
                                             actionButton("profilePlotClear", "Clear Plot!", width = '100%', class = "btn-warning")
                                           ),
                                           mainPanel(
                                             plotOutput("profilePlot")
                                           )
                                  ) # tabPanel "Profile Plots" close
                      )
             ), # tabPanel "Plots" close
             ### OVERVIEW TABLE OVER UP- AND DOWNREGULATED GENES ### 
             tabPanel("Differential Expression",
                      sidebarPanel(
                        h4("Specify Contrast"),
                        selectInput("contrastUpDown_1", "First Experimental Group:", "-"),
                        selectInput("contrastUpDown_2", "Second Experimental Group:", "-"),
                        textOutput("foldChangeInfo"),
                        actionButton("addToOverview", "Add to table!", width = '100%',class = "btn-warning"),
                        actionButton("clearOverview", "Clear table!", width = '100%',class = "btn-warning"),
                      ), # sidebarPanel close
                      mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel("Table", textOutput("overviewInfo"), downloadButton("downloadOverview", "Download"), dataTableOutput("overviewTable"), 
                                             bsModal("diffExpressionResults", "Differential Expression Results", "will_be_triggered_manually", size = "large",
                                                     tabsetPanel(type = "tabs",
                                                                 tabPanel("Gene Expression", tabsetPanel(type = "tabs",
                                                                                                 tabPanel("All Genes", downloadButton("downloadAllResults", "Download"), dataTableOutput("diffResults")),
                                                                                                 tabPanel("Significant Genes", downloadButton("downloadSignResults", "Download"), dataTableOutput("signDiffResults"))
                                                                                                )
                                                                          ),
                                                                 tabPanel("Plots", tabsetPanel(type = "tabs",
                                                                                               #tabPanel("Heatmap", plotOutput("tbc")),
                                                                                               tabPanel("Volcano Plot",
                                                                                                        sidebarPanel(
                                                                                                          h4("Select LogFC Threshold:"),
                                                                                                          sliderInput("volcanoFcThreshold",
                                                                                                                      "Threshold (absolute)",
                                                                                                                      min = 0,
                                                                                                                      max = 3,
                                                                                                                      step = 0.1,
                                                                                                                      value = 1
                                                                                                                      ),
                                                                                                          actionButton("volcPlotButton", "Refresh Plot!", width = '100%', class = "btn-warning")
                                                                                                        ),
                                                                                                        mainPanel(plotOutput("volcanoPlot", brush = "volcanoBrush"), verbatimTextOutput("volcanoInfo"))
                                                                                                        )
                                                                                               )
                                                                          )
                                                                 ) # Modal tabsetPanel close
                                                     ) # Modal close
                                             ),
                                    tabPanel("Venn Diagram", plotOutput("vennDiagram")),
                                    tabPanel("UpSet Plot", plotOutput("upsetPlot"))
                                    
                          ) # tabsetPanel of mainPanel close
                        ) # mainPanel close
             ) # tabPanel "Overview" close
  ) # navBarPage close
) # ui close