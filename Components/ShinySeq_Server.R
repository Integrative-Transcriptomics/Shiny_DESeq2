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
    withProgress(message = "Reading Data", detail = "Reading raw counts", value = 0, {
      countFile = input$countFile
      rawdat = read.table(countFile$datapath, header = TRUE)
      
      ## INPUT INFO DATA ## 
      incProgress(0.2, detail = "Reading experimental design")
      infoFile = input$infoFile
      infodat = read.csv(infoFile$datapath, sep = '\t')
      
      ## INPUT GFF FILE ## 
      incProgress(0.2, detail = "Reading GFF")
      gffFilePath = input$gffFile
      gffdat <<- checkGFF(as.data.frame(readGFF(gffFilePath$datapath)))
      
      
      ## SORT DATA ## 
      incProgress(0.2, detail = "Parsing and disyplaing data")
      dats = sortThatData(rawdat, infodat, gffdat)
      rawCounts <<- dats[[1]]
      infoData <<- dats[[2]]
      
      ## DISPLAY DATA ##
      output$countTable = renderDataTable(rawCounts, rownames = TRUE)
      output$designTable = renderDataTable(infoData, rownames = TRUE)
      incProgress(0.4, detail = "Done.")
      
      ## UPDATE SELECINPUT FOR DESEQ BASED ON INFO DATA ##
      updateSelectInput(session, "variable", choices = factor(colnames(infoData)))
    })
  }) # upload button close
  
  # download counts
  output$downloadCounts = downloadHandler(
    filename = function() {
      paste0("sorted_counts.txt")
    },
    content = function(file) {
      write.csv(rawCounts, file, row.names = TRUE)
    }
  )
    
  # download design
  output$downloadDesign = downloadHandler(
    filename = function() {
      paste0("sorted_design.txt")
    },
    content = function(file) {
      write.csv(infoData, file, row.names = TRUE)
    }
  )
    
  
  ##################
  ### RUN DESEQ2 ###
  ##################
  observeEvent(input$analyze, {
    req(input$upload)

    if(length(levels(as.factor(infoData[,input$variable]))) == nrow(infoData)){
        # selecting a column that contains no replicates results in crash of DESeq
        showNotification("Error: The design matrix has the same number of samples and coefficients to fit,
                          so estimation of dispersion is not possible. Please select a column that contains replicates.", type = "error")
    }
    else{
      withProgress(message = "Running DESeq", detail = "Creating DESeq-Dataset", value = 0, {
        ## CREATE DESEQ DATASET ##
        designFormula = as.formula(paste0("~",input$variable))
        incProgress(0.2, detail = "Creating DESeq-Dataset")
        dds <<- DESeqDataSetFromMatrix(countData = rawCounts[-1], colData =infoData, design = designFormula) 
        dds <<- DESeq(dds)                                    
        ## DISPLAY NORMALIZED COUNTS ##
        incProgress(0.2, detail = "Applying specified normalization")
        # normalization based on selected method (radio button)
        if(input$normMethod == "Size Factor Division"){
          normCounts <<- counts(dds, normalized = TRUE)
        }
        else if(input$normMethod == "VST"){
          normObject <<- varianceStabilizingTransformation(dds)  # S4 object, will e.g. be used in PCA
          normCounts <<- assay(normObject)
        }
        # else if(input$normMethod == "Quantile Normalization"){
        #   # Quantile normalization is not included in DESeq2 => log-transform counts => remove -inf-values => normalize
        #   logCount = logTransform(dats[[1]][,-1])
        #   normCounts = normalize.quantiles(as.matrix(logCount), copy = TRUE)
        # }
          
        # TPM-normalization (always performed)
        incProgress(0.2, detail = "Applying TPM normalization")
        tpmTable <<- normalizeTPM(rawCounts, gffdat)
        
        incProgress(0.2, detail = "Rendering normalized tables")
        output$tpmTable = renderDataTable(datatable(tpmTable) %>% formatRound(columns = c(1:ncol(tpmTable)), digits = 2), rownames = TRUE)
        output$normalizedTable = renderDataTable(datatable(normCounts) %>% formatRound(columns = c(1:ncol(normCounts)), digits = 2), rownames = TRUE)
        
        incProgress(0.2, detail = "Done.")
      })
      # Log-transform (if required). Will be used for heatmaps
      log2normCounts <<- normCounts
      # if data was normalized by size factor division => data is not yet log transformed!
      if(input$normMethod == "Size Factor Division"){
        log2normCounts <<- logTransform(normCounts)
      }
      
      # Update PCA select inputs: 
      is_treatment = grepl("condition", colnames(colData(dds)), ignore.case = TRUE)
      updateSelectInput(session, "pca1", choices = colnames(colData(dds))[is_treatment])
      updateSelectInput(session, "pca2", choices = c("None", colnames(colData(dds))[is_treatment]))
      
      # Update select inputs for differential expression section:
      variables <<- levels(factor(infoData[, c(input$variable)]))
      #updateSelectInput(session, "contrastUpDown_1", choices = variables)
      updateSelectizeInput(session, "contrastUpDown_1", choices = variables)
      
      # Update selectize input for profile plots (gene names):
      updateSelectizeInput(session, "profileGenes", choices = gsub(".*, ", "", row.names(tpmTable)))
      
    }
  }) # analyze button close
  
  # Download normalized counts
  output$downloadNormalizedCounts = downloadHandler(
    filename = function() {
      paste0("normalized_counts.txt")
    },
    content = function(file) {
      write.csv(normCounts, file, row.names = TRUE)
    }
  )
        
  # Download TPM counts
  output$downloadTPMCounts = downloadHandler(
    filename = function() {
      paste0("TPM_counts.txt")
    },
    content = function(file) {
      write.csv(tpmTable, file, row.names = TRUE)
    }
  )
        

  #############
  #### PCA ####
  #############
  
  observeEvent(input$pcaPlot, {
    req(input$pca1)
          
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
    output$pca = renderPlot({makePCA(pcaData = pca, pcaGroups = pcaGroups, fontSize = input$pcaFont)},
      height = input$pcaHeight, width = input$pcaWidth
    ) # render pca plot close
    
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
    req(input$pca1) # using the updating process of pca1 to check if user uploaded and normalized data before (otherwise there would be a crash)
    
    samp_dist = dist(t(log2normCounts))
    #color_gradient = colorRampPalette(c("white", "yellow", "orange" ,"red"))(1000)
    color_gradient = colorRampPalette(c("white", "yellow","orange", "red", "darkred"))(1000)
    plot_experiments = pheatmap(as.matrix(samp_dist), color = color_gradient, silent = TRUE, fontsize = input$experimentHeatFont)
    output$heatExp = renderPlot({plot_experiments}, height = input$experimentHeatHeight, width = input$experimentHeatWidth)
  })
        
        
  # ===== Heatmap of Genes based on highest variance ===== #
        
  observeEvent(input$plotGeneHeat, {
    req(input$pca1)
    
    # Selection of genes:
    numberOfGenes = input$geneHeatNo                                                       # number of genes the user wants to display
    highVarIndex = head(order(rowVars(log2normCounts), decreasing = TRUE), numberOfGenes)  # indexes of the [numberOfGenes] with highest variance
    topVarGenes = log2normCounts[highVarIndex, ]                                           
    topVarGenes = topVarGenes - rowMeans(topVarGenes)                                      
    # column annotation:
    colAnno = data.frame(colData(dds)[, c(input$variable)])                                # annotation is based on the experimental variables the user chose before pressing the analyze button!
    row.names(colAnno) = row.names(colData(dds))                                           # if only one variable is selected, R omits the rownames meaning they need to be re-specified!
    colnames(colAnno) = input$variable                                                     # format column name
    plotGenes = pheatmap(topVarGenes, annotation_col = colAnno, silent = TRUE, annotation_names_col = FALSE, fontsize = input$geneHeatFont, treeheight_row = 5*input$geneHeatFont)
    # plot heatmap:
    output$heatGene = renderPlot({plotGenes}, height = input$geneHeatHeight, width = input$geneHeatWidth)
  })
  
  ###################
  ## PROFILE PLOTS ##
  ###################
  
  observeEvent(input$profilePlotButton, {
    req(input$pca1)
    if(length(input$profileGenes) < 1){
      showNotification("No genes selected for profile plot.", type = "warning")
    }
    else{
      # Function performs quantile normalization of log(tpmTable):
      profilePlots <<- makeProfilePlots(tpmTable, 
                                      geneList = input$profileGenes, 
                                      summarize.replicates = input$averageProfileReplicates,
                                      errorbars = input$profileErrorbars,
                                      fontsize = input$profileFont)  
      
      output$profilePlotColored = renderPlot(profilePlots[[1]], height = input$profileHeight, width = input$profileWidth)
      output$profilePlotBlack = renderPlot(profilePlots[[2]], height = input$profileHeight, width = input$profileWidth)
      
      # Error bars are only possible for averaged sample. Let user know:
      if(!input$averageProfileReplicates & input$profileErrorbars){
        showNotification("NOTE: Error bars are only possible for averaged replicates.", type = "warning")
      }
    }
  })
  
  observeEvent(input$profilePlotClear, {
    req(input$pca1)
    updateSelectizeInput(session, "profileGenes", choices = gsub(".*, ", "", row.names(tpmTable)))
  })
  
  # Download gene profile data
  output$downloadProfileData = downloadHandler(
    filename = function() {
      paste0("profile_data.csv")
    },
    content = function(file) {
      write.csv(profilePlots[[1]][["data"]], file, row.names = TRUE)
    }
  )
  
  
        
  #####################################
  ## DIFFERENTIAL EXPRESSION SECTION ##
  #####################################
  
  observeEvent(input$contrastUpDown_1, {
    if(exists("variables")){
      # Update the 2nd select input based on elements that are not already chosen in the 1st select input
      #updateSelectInput(session, "contrastUpDown_2", choices = c(variables[variables != input$contrastUpDown_1], "Add all"))
      updateSelectInput(session, "contrastUpDown_2", choices = variables[!(variables %in% input$contrastUpDown_1)])
    }
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
  output$plotSizeInfo = renderText("Note: Plot size is updated when entries are added to table.")
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
    else{
      # set 2nd variable as vector so it can be used in loop
      # if(input$contrastUpDown_2 == "Add all"){
      #   secondVariable = variables[variables != input$contrastUpDown_1]
      #   }
      # else{
      #   secondVariable = input$contrastUpDown_2
      # }
      withProgress(message = "Adding new element(s) to table", detail = paste("1 of", length(input$contrastUpDown_1)),value = 0, {
        counter = 1
        for(i in input$contrastUpDown_1){
          counter = counter + 1
          # DESeq crashes if experimental groups are the same
          if(i == input$contrastUpDown_2){
            showNotification(paste(i, "VS", input$contrastUpDown_2, "was ommited, because the experimental groups must differ!"), type = "warning")
            incProgress(1/length(input$contrastUpDown_1), detail = paste(counter, "of", length(input$contrastUpDown_1)))
            next
          }
          # Check for duplicated entries (and ignore them)
          if(paste(i, "VS", input$contrastUpDown_2) %in% overview$data$'Conditions/Comparison'){
            showNotification(paste(i, "VS", input$contrastUpDown_2, "was already added to table and is therefore omitted."), type = "warning")
            incProgress(1/length(input$contrastUpDown_1), detail = paste(counter, "of", length(input$contrastUpDown_1)))
          }
          else{
            # get results, add gene names, product description, fold change, avgTPM and sort
            resultsTable = as.data.frame(results(dds, alpha = input$alpha, contrast = c(input$variable, i, input$contrastUpDown_2)))
            resultsTable = extendAndSortResults(resultsData = resultsTable, gffFile = gffdat, tpmData = tpmTable, contrast1 = i, contrast2 = input$contrastUpDown_2)
            # add to list: 
            resultsList[[length(resultsList)+1]] <<- resultsTable
            #  filter out non-significant (p > alpha, log2FC < 1), get overview (amount of up-/downregulated genes)
            significant_results = filterSignificantGenes(dds_results = resultsTable, alpha = input$alpha, logFCThreshold = 1)
            signResultsList[[length(signResultsList)+1]] <<- significant_results
            significant_overview = significantOverview(significant_results, i, input$contrastUpDown_2)
            # Update overview and genelist, render table
            updateOverviewDataTable(overviewReactiveValues = overview, significantOverviewEntry = significant_overview)
            geneList[[length(geneList)+1]] <<- row.names(significant_results)
            incProgress(1/length(input$contrastUpDown_1), detail = paste(counter, "of", length(input$contrastUpDown_1)))
          }
        } # for-loop-close
        setProgress(value = 1, message = "Adding new element(s) to table", detail = "Done.")
      }) # Progress bar close

      # render plots (upset & venn):
      if(length(geneList) >= 2){
         if(table(overview$data$TOTAL == 0)[1] >= 2){
           # UpsetR will crash if there are are less than two non-empty elements in list
          output$upsetPlot = renderPlot({makeUpset(geneList, overview$data, fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
        }
              
        if(length(geneList) <= 4){
           # Venn diagram with >4 dimensions is too confusing
           output$vennDiagram = renderPlot({makeVenn(geneList, overview$data, fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
        }
         else{
          # if there are >4 entries, make sure the first 4 entries are contained in diagram
          output$vennDiagram = renderPlot({makeVenn(geneList[1:4], overview$data[1:4,], fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
         }
       }
      if(length(geneList) > 4){
        # inform user when maximum 
        showNotification("Warning: Venn diagram only supports 2-4 dimensions. Addition of further dimensions will be ignored.", type = "warning")
      }
    }
  }) # add to overview button close
  

        
  # ===== Differential Expression: Displaying Results (Pop-Up) ===== #
        
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
         # Sort table by absolute logFC before download
         dl_table = as.data.frame(resultsList[[rowIndex]])
         dl_table = dl_table[sort(abs(dl_table$log2FoldChange), decreasing = TRUE, index.return = TRUE)[[2]],]
         write.csv(dl_table, file, row.names = TRUE)
       }
     )
    output$downloadSignResults = downloadHandler(
       filename = function() {
        paste0(overview$data[rowIndex, 2],".csv")
       },
      content = function(file) {
        # Sort table by absolute logFC before download
        dl_table = signResultsList[[rowIndex]]
        dl_table = dl_table[sort(abs(dl_table$log2FoldChange), decreasing = TRUE, index.return = TRUE)[[2]],]
        write.csv(dl_table, file, row.names = TRUE)
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
        output$upsetPlot = renderPlot({makeUpset(geneList, overview$data, fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
        if(length(geneList) <= 4){
          # Venn diagram with >4 dimensions is too confusing
          output$vennDiagram = renderPlot({makeVenn(geneList, overview$data, fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
        }
        else{
           # if there are >4 entries, make sure the first 4 entries are contained in diagram
           output$vennDiagram = renderPlot({makeVenn(geneList[1:4], overview$data[1:4,], fontsize = input$overviewFont)}, height = input$overviewHeight, width = input$overviewWidth)
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
}) # server close
