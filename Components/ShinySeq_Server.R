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
    gffFile = input$gffFile
    gffdat = as.data.frame(readGFF(gffFile$datapath))
    
    ## SORT DATA ## 
    dats = sortThatData(rawdat, infodat, gffdat)
    rawCounts = dats[[1]]
    infoData = dats[[2]]
    
    ## DISPLAY DATA ##
    output$countTable = renderTable(rawCounts, rownames = TRUE)
    output$designTable = renderTable(infoData, rownames = TRUE)
    
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
      
      # else if(input$normMethod == "Quantile Normalization"){
      #   # Quantile normalization is not included in DESeq2 => log-transform counts => remove -inf-values => normalize
      #   logCount = logTransform(dats[[1]][,-1])
      #   normCounts = normalize.quantiles(as.matrix(logCount), copy = TRUE)
      # }
      
      output$normalizedTable = renderTable(normCounts, rownames = TRUE) 
      
      
      ## DISPLAY RESULTS ##
      
      ## Update results select inputs
      variables = levels(factor(infoData[, c(input$variable)]))
      updateSelectInput(session, "contrast1", choices = variables)
      updateSelectInput(session, "contrast2", choices = variables)
      
      observeEvent(input$results, {
        if(is.null(dds)){
          showNotification("Please run DESeq first", type = "error")
        }
        req(dds)
        
        # get results
        ddsRes = results(dds, alpha = input$alpha, contrast = c(input$variable, input$contrast1, input$contrast2)) # significance level is chosen by user via slider
        output$resTable = renderTable(as.data.frame(ddsRes), rownames = TRUE)
        
        significant_results = filterSignificantGenes(dds_results = ddsRes, alpha = input$alpha, logFCThreshold = 1)
        output$significantTable = renderTable(as.data.frame(significant_results), rownames = TRUE)
        
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
        updateSelectInput(session, "pca2", choices = colnames(colData(dds))[is_treatment])
        
        observeEvent(input$pcaPlot, {
          # plot PCA based on chosen normalization method
          if(input$normMethod == "Size Factor Division"){
            pca = plotPCA(rlog(dds), intgroup = c(input$pca1, input$pca2)) 
          }
          else{
            # get data from plotPCA so it can be modified
            pca = plotPCA(normObject, intgroup = c(input$pca1, input$pca2))
          }
          # plot
          output$pca = renderPlot({
            ggplot(pca[["data"]], aes(x = PC1, y = PC2, color = pca[["data"]][[input$pca1]], shape = pca[["data"]][[input$pca2]])) +
              geom_point(size = 2) +
              theme(legend.title = element_blank()) + 
              labs(x = pca[["labels"]][["x"]], y = pca[["labels"]]["y"])
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
