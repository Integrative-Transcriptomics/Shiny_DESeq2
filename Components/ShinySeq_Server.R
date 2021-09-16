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
    output$countTable = renderTable(rawCounts, rownames = TRUE)
    output$designTable = renderTable(infoData, rownames = TRUE)
    
    # download counts
    output$downloadCounts = downloadHandler(
      filename = function() {
        paste0("sorted_counts.txt")
      },
      content = function(file) {
        write.csv(dats[[1]], file, row.names = TRUE)
      }
    )
    
    ## UPDATE SELECINPUT FOR DESEQ BASED ON INFO DATA ##
    updateSelectInput(session, "variable", choices = factor(colnames(infoData)))
    
    ##################
    ### RUN DESEQ2 ###
    ##################
    observeEvent(input$analyze, {

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
      # TPM-normalization (always performed)
      tpmTable = normalizeTPM(rawCounts, gffdat)
      output$tpmTable = renderTable(tpmTable, rownames = TRUE)
      
      # else if(input$normMethod == "Quantile Normalization"){
      #   # Quantile normalization is not included in DESeq2 => log-transform counts => remove -inf-values => normalize
      #   logCount = logTransform(dats[[1]][,-1])
      #   normCounts = normalize.quantiles(as.matrix(logCount), copy = TRUE)
      # }
      
      output$normalizedTable = renderTable(normCounts, rownames = TRUE) 
      
      ## Update results select inputs
      variables = levels(factor(infoData[, c(input$variable)]))
      updateSelectInput(session, "contrast1", choices = variables)
      updateSelectInput(session, "contrast2", choices = variables)
      
      updateSelectInput(session, "contrastUpDown_1", choices = variables)
      updateSelectInput(session, "contrastUpDown_2", choices = variables)
      
      #########################
      ## UP-/DOWN-REGULATION ##
      #########################
      
      # Default Setting:
      overview = reactiveValues(data = data.frame("Conditions/Comparison" = character(0), 
                                                  "UP" = numeric(0),
                                                  "DOWN" = numeric(0),
                                                  "TOTAL" = numeric(0),
                                                  "Actions" = shinyInput(actionButton, 0, 'button_', label = "Delete", onclick = 'Shiny.onInputChange(\"delete_button\",  this.id)')
                                                  )
                                )  # empty table, will get updated
      output$overviewTable = renderDataTable(overview$data, escape = FALSE)
      output$overviewInfo = renderText("Add at least 2 sets to overview table in order to display venn diagram and upset plot")
      geneList = list()  # for venn diagram and UpSet plot
      output$vennDiagram = renderPlot({ggplot()})
      output$upsetPlot = renderPlot({ggplot()})
      
      # Adding things to overview table (as well as Venn & UpSet Plot)
      observeEvent(input$addToOverview, {
        if(is.null(dds)){
          showNotification("Please run DESeq first", type = "error")
        }
        req(dds)
        
        if(input$contrastUpDown_1 == input$contrastUpDown_2){
          # DESeq crashes if experimental groups are the same
          showNotification("Experimental groups must differ!", type = "error")
        }
        else if(paste(input$contrastUpDown_1, "VS", input$contrastUpDown_2) %in% overview$data[,1]){
          showNotification("Contrast was already added to overview table!", type = "error")
        }
        else{
          # get results, filter out non-significant (p > alpha, log2FC < 1), get overview (amount of up-/downregulated genes)
          resultsOverview = results(dds, alpha = input$alpha, contrast = c(input$variable, input$contrastUpDown_1, input$contrastUpDown_2))
          significant_results = filterSignificantGenes(dds_results = resultsOverview, alpha = input$alpha, logFCThreshold = 1)
          significant_overview = significantOverview(significant_results, input$contrastUpDown_1, input$contrastUpDown_2)
          
          # Update & render table
          overview$data = rbind(overview$data[,-5], significant_overview)
          overview$data$Actions = shinyInput(actionButton, nrow(overview$data), 'button_', label = "Delete", onclick = 'Shiny.onInputChange(\"delete_button\",  this.id)')
          
          # Update list & render venn Diagram and UpSet plot
          geneList[[length(geneList)+1]] <<- row.names(significant_results)
          if(length(geneList) >= 2){
            output$upsetPlot = renderPlot({makeUpset(geneList, overview$data)})
            if(length(geneList) <= 7){
              # ggVennDiagram only supports 2-7 dimensions -> ignore >7 dimensions
              output$vennDiagram = renderPlot({makeVenn(geneList, overview$data)})
            }
          }
          if(length(geneList) > 7){
            # inform user when maximum 
            showNotification("Warning: Venn diagram only supports 2-7 dimensions. Addition of further dimensions will be ignored.", type = "warning")
            }
          }
      }
      ) # add contrast to overview table close
      
      # Download button overview table as .tsv:
      output$downloadOverview = downloadHandler(
        filename = "overview.tsv",
        content = function(file) {
          write.table(overview$data[,-5], file, row.names = FALSE, sep = "\t")
        }
      )
      
      # Clear specific row:
      observeEvent(input$delete_button, {
        print(paste("pressed button", input$delete_button))
        overview$data <<- overview$data[,-5]
        # Update overview table
        rowIndex = as.numeric(strsplit(input$delete_button, "_")[[1]][2])
        overview$data <<- overview$data[-rowIndex,]  # also removes column with 'delete'-buttons so a new column with updated IDs of action buttons can be added
        overview$data$Actions <<- shinyInput(actionButton, nrow(overview$data), 'button_', label = "Delete", onclick = 'Shiny.onInputChange(\"delete_button\",  this.id)')
        # Update venn diagram & Upset: 
        geneList[[rowIndex]] <<- NULL
        if(length(geneList) >= 2){
          output$upsetPlot = renderPlot({makeUpset(geneList, overview$data)})
          if(length(geneList) <= 7){
            # ggVennDiagram only supports 2-7 dimensions -> ignore >7 dimensions
            output$vennDiagram = renderPlot({makeVenn(geneList, overview$data)})
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
        overview$data = data.frame("Conditions/Comparison" = character(0), 
                                   "UP" = numeric(0),
                                   "DOWN" = numeric(0),
                                   "TOTAL" = numeric(0),
                                   "Actions" = shinyInput(actionButton, 0, 'button_', label = "Delete", onclick = 'Shiny.onInputChange(\"delete_button\",  this.id)')
                        )
        # clear list for venn and UpSet
        geneList <<- list()
        output$vennDiagram = renderPlot({ggplot()})
        output$upsetPlot = renderPlot({ggplot()})
      })
      
      
      ## DISPLAY RESULTS ##
      observeEvent(input$results, {
        if(is.null(dds)){
          showNotification("Please run DESeq first", type = "error")
        }
        req(dds)
        
        # get results, add gene-, description-, foldchange and average-TPM (log and absolute) -columns
        ddsRes = as.data.frame(results(dds, alpha = input$alpha, contrast = c(input$variable, input$contrast1, input$contrast2))) # significance level is chosen by user via slider
        ddsRes = addGeneNameCol(ddsRes)
        ddsRes = addDescriptionCol(ddsRes, gffdat)  # read-in of gff-data already happened at the very beginning
        ddsRes = addFoldChangeCol(ddsRes)
        ddsRes = addAverageTPM(ddsRes, tpmTable)
        # sort and render table: 
        ddsRes = ddsRes[,c(7,1,2,9,3,10,4:6,8)]  # sort
        output$resTable = renderTable(ddsRes, rownames = TRUE)
        
        significant_results = filterSignificantGenes(dds_results = ddsRes, alpha = input$alpha, logFCThreshold = 1)
        output$significantTable = renderTable(significant_results, rownames = TRUE)
        
        # in case user specified a different normalization method, display message that results table and volcano plots are still based on DESeq standard:
        if(input$normMethod != "Size Factor Division"){
          output$resText_all = renderText("NOTE: Differential Expression Results will always be based on DESeq2's standard normalization method (Size Factor Division).")
          output$resText_sig = renderText("NOTE: Differential Expression Results will always be based on DESeq2's standard normalization method (Size Factor Division).")
        }
        
        # download results
        output$downloadResults = downloadHandler(
          filename = function() {
            paste0(input$contrast1,"_vs_",input$contrast2,".csv")
          },
          content = function(file) {
            write.csv(ddsRes, file, row.names = TRUE)
          }
        )
        #download significant results
        output$downloadSignificant = downloadHandler(
          filename = function() {
            paste0(input$contrast1,"_vs_",input$contrast2,".csv")
          },
          content = function(file) {
            write.csv(significant_results, file, row.names = TRUE)
          }
        )
        
        
        ##################
        ##### PLOTS ######
        ##################
        
        #### PCA ####
        # update PCA select inputs: 
        is_treatment = grepl("condition", colnames(colData(dds)), ignore.case = TRUE)
        updateSelectInput(session, "pca1", choices = colnames(colData(dds))[is_treatment])
        updateSelectInput(session, "pca2", choices = c(colnames(colData(dds))[is_treatment], "-"))
        
        
        observeEvent(input$pcaPlot, {
          
          # set variables for PCA:
          if(input$pca2 == "-"){
            pcaGroups = input$pca1
          }
          else{
            pcaGroups = c(input$pca1, input$pca2)
          }
          
          # plot PCA based on chosen normalization method
          if(input$normMethod == "Size Factor Division"){
            pca = plotPCA(rlog(dds), intgroup = pcaGroups) 
          }
          else{
            # get data from plotPCA so it can be modified
            pca = plotPCA(normObject, intgroup = pcaGroups)
          }
          # plot
          output$pca = renderPlot({
            # ggplot(pca[["data"]], aes(x = PC1, y = PC2, color = pca[["data"]][[input$pca1]], shape = pca[["data"]][[input$pca2]])) +
            #   geom_point(size = 2) +
            #   theme(legend.title = element_blank()) + 
            #   labs(x = pca[["labels"]][["x"]], y = pca[["labels"]]["y"])
            makePCA(pcaData = pca, pcaGroups = pcaGroups)
          }) # render pca plot close
          
          # pca interactive brush info
          output$pca_info = renderPrint({
            brushedPoints(pca[["data"]], input$pca_brush)
          })
        }) # PCA button close
        
        
        #### VOLCANO PLOT ####
        observeEvent(input$volcPlot, {
          # in case user specified a different normalization method, display message that results table and volcano plots are still based on DESeq standard:
          if(input$normMethod != "Size Factor Division"){
            showNotification("NOTE: Volcano plot is based on Size Factor Divison.")
          }
          volc_plot_and_data = erupt(ddsRes, input$volcFcThr, input$alpha)  # creates plot object ([[1]]) and data ([[2]])
          
          # render plot:
          output$volc = renderPlot({volc_plot_and_data[[1]]}) 
          # volcano plot interactive brush info:
          output$volc_info = renderPrint({
            brushedPoints(volc_plot_and_data[[2]], input$volc_brush)
          })
        }) # volcano plot button close
        
        
        #### HEATMAP OF EXPERIMENTS ####
        log2normCounts = normCounts
        # if data was normalized by size factor division => data is not yet log transformed!
        if(input$normMethod == "Size Factor Division"){
          log2normCounts = logTransform(normCounts)
        }
        # render plot
        samp_dist = dist(t(log2normCounts))
        #color_gradient = colorRampPalette(c("white", "yellow", "orange" ,"red"))(1000)
        color_gradient = colorRampPalette(c("white", "yellow","orange", "red", "darkred"))(1000)
        plot_experiments = pheatmap(as.matrix(samp_dist), color = color_gradient)
        output$heatExp = renderPlot({plot_experiments})
        
        #### HEATMAP OF TOP GENES (BASED ON VARIANCE) ####
        observeEvent(input$plotGeneHeat, {
          # Selection of genes:
          numberOfGenes = input$geneHeatNo                                                       # number of genes the user wants to display
          highVarIndex = head(order(rowVars(log2normCounts), decreasing = TRUE), numberOfGenes)  # indexes of the [numberOfGenes] with highest variance
          topVarGenes = log2normCounts[highVarIndex, ]                                           # subset log2 transformed dataset accordingly
          topVarGenes = topVarGenes - rowMeans(topVarGenes)                                      # mean centering to acquire (log2-)deviation from the mean
          # column annotation:
          colAnno = data.frame(colData(dds)[, c(input$variable)])                                # annotation is based on the experimental variables the user chose before pressing the analyze button!
          row.names(colAnno) = row.names(colData(dds))                                           # if only one variable is selected, R omits the rownames meaning the need to be re-specified!
          colnames(colAnno) = input$variable                                                     # format column name
          plotGenes = pheatmap(topVarGenes, annotation_col = colAnno)
          # plot heatmap:
          output$heatGene = renderPlot({plotGenes}, height = input$geneHeatHeight)
          
        }) # gene Heatmap button close
        
        #### BOXPLOTS ####
        output$boxplot = renderPlot(boxplot(log2normCounts))
      }) # results button close
    }) # analyze button close
  }) # upload button close
}) # server close
