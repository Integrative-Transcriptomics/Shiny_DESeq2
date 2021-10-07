#### ========= SHINYSEQ SERVER ========== ####
server = shinyServer(function(input, output, session){
  options(shiny.sanitize.errors = TRUE)
  
  observeEvent(input$analyze, {
    
    ## Error message if data has not been uploaded before ## 
    if(is.null(input$countFile) | is.null(input$infoFile) | is.null(input$gffFile)){
      showNotification("Please upload count data, sample preparation info and .gff-file", type = "error")
    }
  })
  
  
  ##############
  ### UPLOAD ###
  ##############
  ## Display Data if user presses 'Upload!' button ## 
  observeEvent(input$upload,{
    
    ## Exception handling and display of error message if data is missing ##
    if(is.null(input$countFile) | is.null(input$infoFile) | is.null(input$gffFile)){
      showNotification("Please upload count data, sample preparation info and .gff-file", type = "error")
    }
    req(input$countFile)
    req(input$infoFile)
    req(input$gffFile)
    
    ## INPUT RAW DATA ##
    countFile = input$countFile
    rawdat = read.table(countFile$datapath, header = TRUE)
    
    ## INPUT INFO DATA ## 
    infoFile = input$infoFile
    infodat = read.csv(infoFile$datapath, sep = '\t')
    
    ## INPUT GFF FILE ## 
    gffFilePath = input$gffFile
    gffdat = as.data.frame(readGFF(gffFilePath$datapath))
    
    ## SORT DATA ## 
    dats = sortThatData(rawdat, infodat, gffdat)
    rawCounts = dats[[1]]
    infoData = dats[[2]]
    
    ## DISPLAY DATA ##
    output$countTable = renderDataTable(rawCounts, rownames = TRUE)
    output$designTable = renderDataTable(infoData, rownames = TRUE)
    
    # download counts
    output$downloadCounts = downloadHandler(
      filename = function() {
        paste0("sorted_counts.txt")
      },
      content = function(file) {
        write.csv(dats[[1]], file, row.names = TRUE)
      }
    )
    
    # download design
    output$downloadDesign = downloadHandler(
      filename = function() {
        paste0("sorted_design.txt")
      },
      content = function(file) {
        write.csv(dats[[2]], file, row.names = TRUE)
      }
    )
    
    ## UPDATE SELECINPUT FOR DESEQ BASED ON INFO DATA ##
    updateSelectInput(session, "variable", choices = factor(colnames(infoData)))
    
    ##################
    ### RUN DESEQ2 ###
    ##################
    observeEvent(input$analyze, {
      if(length(levels(as.factor(infoData[,input$variable]))) == nrow(infoData)){
        # selecting a column that contains no replicates results in crash of DESeq
        showNotification("Error: The design matrix has the same number of samples and coefficients to fit,
                          so estimation of dispersion is not possible. Please select a column that contains replicates.", type = "error")
      }
      else{
        ## CREATE DESEQ DATASET ##
        designFormula = as.formula(paste0("~",input$variable))
        dds = DESeqDataSetFromMatrix(countData = rawCounts[-1], colData =infoData, design = designFormula) 
        dds = DESeq(dds)                                    
        
        ## DISPLAY NORMALIZED COUNTS ##
        # normalization based on selected method (radio button)
        if(input$normMethod == "Size Factor Division"){
          normCounts = counts(dds, normalized = TRUE)
        }
        else if(input$normMethod == "VST"){
          normObject = varianceStabilizingTransformation(dds)  # S4 object, will e.g. be used in PCA
          normCounts = assay(normObject)
        }
        # else if(input$normMethod == "Quantile Normalization"){
        #   # Quantile normalization is not included in DESeq2 => log-transform counts => remove -inf-values => normalize
        #   logCount = logTransform(dats[[1]][,-1])
        #   normCounts = normalize.quantiles(as.matrix(logCount), copy = TRUE)
        # }
        output$normalizedTable = renderDataTable(datatable(normCounts) %>% formatRound(columns = c(1:ncol(normCounts)), digits = 2), rownames = TRUE)
        
        # Download normalized counts
        output$downloadNormalizedCounts = downloadHandler(
          filename = function() {
            paste0("normalized_counts.txt")
          },
          content = function(file) {
            write.csv(normCounts, file, row.names = TRUE)
          }
        )
        
        # TPM-normalization (always performed)
        tpmTable = normalizeTPM(rawCounts, gffdat)
        output$tpmTable = renderDataTable(datatable(tpmTable) %>% formatRound(columns = c(1:ncol(tpmTable)), digits = 2), rownames = TRUE)
        
        # Download TPM counts
        output$downloadTPMCounts = downloadHandler(
          filename = function() {
            paste0("TPM_counts.txt")
          },
          content = function(file) {
            write.csv(tpmTable, file, row.names = TRUE)
          }
        )
        
        
        # Log-transform (if required). Will be used for heatmaps
        log2normCounts = normCounts
        # if data was normalized by size factor division => data is not yet log transformed!
        if(input$normMethod == "Size Factor Division"){
          log2normCounts = logTransform(normCounts)
        }
        
        #############
        #### PCA ####
        #############
        
        # update PCA select inputs: 
        is_treatment = grepl("condition", colnames(colData(dds)), ignore.case = TRUE)
        updateSelectInput(session, "pca1", choices = colnames(colData(dds))[is_treatment])
        updateSelectInput(session, "pca2", choices = c("None", colnames(colData(dds))[is_treatment]))
        
        observeEvent(input$pcaPlot, {
          
          # set variables for PCA:
          if(input$pca2 == "None"){
            pcaGroups = input$pca1
          }
          else{
            pcaGroups = c(input$pca1, input$pca2)
          }
          
          # plot PCA based on chosen normalization method
          if(input$normMethod == "Size Factor Division"){
            # get data from plotPCA so it can be modified
            pca = plotPCA(rlog(dds), intgroup = pcaGroups) 
          }
          else{
            # get data from plotPCA so it can be modified
            pca = plotPCA(normObject, intgroup = pcaGroups)
          }
          # plot
          output$pca = renderPlot({
            makePCA(pcaData = pca, pcaGroups = pcaGroups)
          }) # render pca plot close
          
          # pca interactive brush info
          output$pca_info = renderPrint({
            brushedPoints(pca[["data"]], input$pca_brush)
          })
        }) # PCA button close
        
        ##################
        #### BOXPLOTS ####
        ##################
        
        output$boxplot = renderPlot(boxplot(log2normCounts))
        
        ################################
        #### HEATMAPS OF EXEPRIMENT ####
        ################################
        
        
        # ===== Heatmap of samples ===== #
        
        observeEvent(input$plotExperimentHeat, {
          samp_dist = dist(t(log2normCounts))
          #color_gradient = colorRampPalette(c("white", "yellow", "orange" ,"red"))(1000)
          color_gradient = colorRampPalette(c("white", "yellow","orange", "red", "darkred"))(1000)
          plot_experiments = pheatmap(as.matrix(samp_dist), color = color_gradient, silent = TRUE)
          output$heatExp = renderPlot({plot_experiments}, height = input$experimentHeatHeight, width = input$experimentHeatWidth)
        })
        
        
        # ===== Heatmap of Genes based on highest variance ===== #
        
        observeEvent(input$plotGeneHeat, {
          # Selection of genes:
          numberOfGenes = input$geneHeatNo                                                       # number of genes the user wants to display
          highVarIndex = head(order(rowVars(log2normCounts), decreasing = TRUE), numberOfGenes)  # indexes of the [numberOfGenes] with highest variance
          topVarGenes = log2normCounts[highVarIndex, ]                                           
          topVarGenes = topVarGenes - rowMeans(topVarGenes)                                      
          # column annotation:
          colAnno = data.frame(colData(dds)[, c(input$variable)])                                # annotation is based on the experimental variables the user chose before pressing the analyze button!
          row.names(colAnno) = row.names(colData(dds))                                           # if only one variable is selected, R omits the rownames meaning they need to be re-specified!
          colnames(colAnno) = input$variable                                                     # format column name
          plotGenes = pheatmap(topVarGenes, annotation_col = colAnno, silent = TRUE)
          # plot heatmap:
          output$heatGene = renderPlot({plotGenes}, height = input$geneHeatHeight)
        })
        
        
        #####################################
        ## DIFFERENTIAL EXPRESSION SECTION ##
        #####################################
        
        # Update select inputs for differential expression
        variables = levels(factor(infoData[, c(input$variable)]))
        updateSelectInput(session, "contrastUpDown_1", choices = variables)
        observeEvent(input$contrastUpDown_1, {
          updateSelectInput(session, "contrastUpDown_2", choices = variables[variables != input$contrastUpDown_1])
        })
        
        
        # ===== Differential Expression: Initialization ===== #
        
        # Default Settings:
        overview = reactiveValues(data = data.frame("Set" = character(0),
                                                    "Conditions/Comparison" = character(0), 
                                                    "UP" = numeric(0),
                                                    "DOWN" = numeric(0),
                                                    "TOTAL" = numeric(0),
                                                    "Genes" = shinyInput(actionButton, 0, 'button_', label = "Show Genes", onclick = 'Shiny.setInputValue(\"genes_button\",  this.id.concat(\"_\", Math.random()))'),
                                                    "Delete" = shinyInput(actionButton, 0, 'button_', label = "Delete", onclick = 'Shiny.setInputValue(\"delete_button\",  this.id.concat(\"_\", Math.random()))')
        )
        )  # empty table, will get updated
        output$overviewTable = renderDataTable(overview$data, rownames = FALSE, escape = FALSE)
        output$foldChangeInfo = renderText("The foldchange is always calculated as first group/second group or log(first group) - log(second group), respectively")
        output$overviewInfo = renderText("Add at least 2 sets to overview table in order to display venn diagram and upset plot")
        output$vennDiagram = renderPlot({ggplot()})
        output$upsetPlot = renderPlot({ggplot()})
        
        # Lists:
        geneList = list()         # for venn diagram and UpSet plot
        resultsList = list()      # will contain all results
        signResultsList = list()  # will contain significant results only
        
        
        # ===== Differential Expression: Adding things to overview table (as well as Venn & UpSet Plot) ===== #
        
        observeEvent(input$addToOverview, {
          if(is.null(dds)){
            showNotification("Please run DESeq first", type = "error")
          }
          req(dds)
          
          if(input$contrastUpDown_1 == input$contrastUpDown_2){
            # DESeq crashes if experimental groups are the same
            showNotification("Experimental groups must differ!", type = "error")
          }
          else if(paste(input$contrastUpDown_1, "VS", input$contrastUpDown_2) %in% overview$data$'Conditions/Comparison'){
            showNotification("Contrast was already added to overview table!", type = "error")
          }
          else{
            # get results, add gene names, product description, fold change, avgTPM and sort
            resultsTable = results(dds, alpha = input$alpha, contrast = c(input$variable, input$contrastUpDown_1, input$contrastUpDown_2))
            resultsTable = extendAndSortResults(resultsData = resultsTable, gffFile = gffdat, tpmData = tpmTable)
            # add to list: 
            resultsList[[length(resultsList)+1]] <<- resultsTable
            #  filter out non-significant (p > alpha, log2FC < 1), get overview (amount of up-/downregulated genes)
            significant_results = filterSignificantGenes(dds_results = resultsTable, alpha = input$alpha, logFCThreshold = 1)
            signResultsList[[length(signResultsList)+1]] <<- significant_results
            significant_overview = significantOverview(significant_results, input$contrastUpDown_1, input$contrastUpDown_2)
            
            # Update & render table
            overview$data = rbind(overview$data[,-c(1,6,7)], significant_overview)
            overview$data$Set = make.unique.3(rep("Set", nrow(overview$data)), sep = "_")
            overview$data = overview$data[,c(5,1:4)]
            overview$data$Genes = shinyInput(actionButton, nrow(overview$data), 'button_', label = "Show Genes", onclick = 'Shiny.setInputValue(\"genes_button\",  this.id.concat(\"_\", Math.random()))')
            overview$data$Delete = shinyInput(actionButton, nrow(overview$data), 'button_', label = "Delete", onclick = 'Shiny.setInputValue(\"delete_button\",  this.id.concat(\"_\", Math.random()))')
            
            # Update list & render venn Diagram and UpSet plot
            geneList[[length(geneList)+1]] <<- row.names(significant_results)
            if(length(geneList) >= 2){
              if(table(overview$data$TOTAL == 0)[1] >= 2){
                # UpsetR will crash if there are are less than two non-empty elements in list
                output$upsetPlot = renderPlot({makeUpset(geneList, overview$data)})
              }
              
              if(length(geneList) <= 4){
                # Venn diagram with >4 dimensions is too confusing
                output$vennDiagram = renderPlot({makeVenn(geneList, overview$data)})
              }
              else{
                # if there are >4 entries, make sure the first 4 entries are contained in diagram
                output$vennDiagram = renderPlot({makeVenn(geneList[1:4], overview$data[1:4,])})
              }
            }
            if(length(geneList) > 4){
              # inform user when maximum 
              showNotification("Warning: Venn diagram only supports 2-4 dimensions. Addition of further dimensions will be ignored.", type = "warning")
            }
          }
        }
        ) # add contrast to overview table close
        
        
        # ===== Differential Expression: Displaying Results ===== #
        
        observeEvent(input$genes_button, {
          rowIndex = as.numeric(strsplit(input$genes_button, "_")[[1]][2])
          # == Tables ==
          output$diffResults = renderDataTable(datatable(as.data.frame(resultsList[[rowIndex]])) %>% formatRound(columns = c(2:9), digits = 2))
          output$signDiffResults = renderDataTable(datatable(as.data.frame(signResultsList[[rowIndex]])) %>% formatRound(columns = c(2:9), digits = 2))
          toggleModal(session, modalId = "diffExpressionResults", toggle = "open")
          
          # == Plots ==
          
          # Heatmap
          
          
          # Volcano
          observeEvent(input$volcPlotButton, {
            # in case user specified a different normalization method, display message that results table and volcano plots are still based on DESeq standard:
            if(input$normMethod != "Size Factor Division"){
              showNotification("NOTE: Volcano plot is based on Size Factor Divison.")
            }
            volc_plot_and_data = erupt(resultsList[[rowIndex]], input$volcanoFcThreshold, input$alpha)  # creates plot object ([[1]]) and data ([[2]])
            
            # render plot:
            output$volcanoPlot = renderPlot({volc_plot_and_data[[1]]}) 
            # volcano plot interactive brush info:
            output$volcanoInfo = renderPrint({
              brushedPoints(volc_plot_and_data[[2]], input$volcanoBrush)
            })
          }) # volcano plot button close
          
          
          # == Table downloads ==
          output$downloadAllResults = downloadHandler(
            filename = function() {
              paste0(overview$data[rowIndex, 2],".csv")
            },
            content = function(file) {
              write.csv(resultsList[[rowIndex]], file, row.names = TRUE)
            }
          )
          output$downloadSignResults = downloadHandler(
            filename = function() {
              paste0(overview$data[rowIndex, 2],".csv")
            },
            content = function(file) {
              write.csv(signResultsList[[rowIndex]], file, row.names = TRUE)
            }
          )
          
        })
        
        
        # ===== Differential Expression: Download Overview Table ===== #
        
        # Overview table as .tsv:
        output$downloadOverview = downloadHandler(
          filename = "overview.tsv",
          content = function(file) {
            write.table(overview$data[,-c(6,7)], file, row.names = FALSE, sep = "\t")
          }
        )
        
        
        # ===== Differential Expression: Deleting things ===== #
        
        # delete (row) - button:
        observeEvent(input$delete_button, {
          # Update overview table
          rowIndex = as.numeric(strsplit(input$delete_button, "_")[[1]][2])
          overview$data <<- overview$data[-rowIndex, -c(6,7)]  # also removes column with 'delete'-buttons so a new column with updated IDs of action buttons can be added
          overview$data$Set = make.unique.3(rep("Set", nrow(overview$data)), sep = "_")
          overview$data$Genes = shinyInput(actionButton, nrow(overview$data), 'button_', label = "Show Genes", onclick = 'Shiny.setInputValue(\"genes_button\",  this.id.concat(\"_\", Math.random()))')
          overview$data$Delete <<- shinyInput(actionButton, nrow(overview$data), 'button_', label = "Delete", onclick = 'Shiny.setInputValue(\"delete_button\",  this.id.concat(\"_\", Math.random()))')
          #row.names(overview$Data) <<- c(1:nrow(overview$Data))
          # Update lists, venn diagram & Upset: 
          geneList[[rowIndex]] <<- NULL
          resultsList[[rowIndex]] <<- NULL
          signResultsList[[rowIndex]] <<- NULL
          if(length(geneList) >= 2){
            output$upsetPlot = renderPlot({makeUpset(geneList, overview$data)})
            if(length(geneList) <= 4){
              # Venn diagram with >4 dimensions is too confusing
              output$vennDiagram = renderPlot({makeVenn(geneList, overview$data)})
            }
            else{
              # if there are >4 entries, make sure the first 4 entries are contained in diagram
              output$vennDiagram = renderPlot({makeVenn(geneList[1:4], overview$data[1:4,])})
            }
          }
          else{
            output$vennDiagram = renderPlot({ggplot()})
            output$upsetPlot = renderPlot({ggplot()})
          }
        })
        
        
        # Clear button for entire overview table:
        observeEvent(input$clearOverview, {
          # clear table:
          #overview$data = overview$data[-(1:nrow(overview$data)),]
          overview$data = data.frame("Set" = character(0),
                                     "Conditions/Comparison" = character(0), 
                                     "UP" = numeric(0),
                                     "DOWN" = numeric(0),
                                     "TOTAL" = numeric(0),
                                     "Genes" = shinyInput(actionButton, 0, 'button_', label = "Show Genes", onclick = 'Shiny.setInputValue(\"genes_button\",  this.id.concat(\"_\", Math.random()))'),
                                     "Delete" = shinyInput(actionButton, 0, 'button_', label = "Delete", onclick = 'Shiny.setInputValue(\"delete_button\",  this.id.concat(\"_\", Math.random()))')
          )
          # clear lists, venn and UpSet
          geneList <<- list()
          resultsList <<- list()
          signResultsList <<- list()
          output$vennDiagram = renderPlot({ggplot()})
          output$upsetPlot = renderPlot({ggplot()})
        })
      }
    }) # analyze button close
  }) # upload button close
}) # server close
