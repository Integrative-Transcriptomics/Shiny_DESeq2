#### =========== SHINYSEQ FUNCTIONS ========== ####
# Functions to be used in the server component

## MISC ##

# Altered version of make.unique
make.unique.2 = function(x, sep='.'){
  # https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
  # Makes unique names for duplicate entries in a vector. Default make.unique leaves the first duplicate unchanged. This function starts enumarting from duplicate 1. 
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}


# Another altered version of make.unique. This version also adds an ID-number when 'x' only contains a single object
make.unique.3 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) >= 1){paste(a, 1:length(a), sep=sep)} else {a}})
}


# Function that enables addition of action buttons to each row of a table
# https://stackoverflow.com/questions/45739303/r-shiny-handle-action-buttons-in-data-table
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}


## RAW DATA MANIPULATION ##

# Many functions require a column 'gbkey' in the GFF. Sometimes the column is called 'type' instead causing several issues.
checkGFF = function(gffFile){
  if(!"gbkey" %in% colnames(gffFile)){
    colnames(gffFile)[colnames(gffFile) == 'type'] <- 'gbkey'
  }
  return(gffFile)
}


# Helper-function that is called in sortThatData(). Concatenates every possible combination of columns containing 'condition' 
addInteractionColumns = function(conditionsColumns, designTable){
  colNumber = ncol(conditionsColumns)
  if(colNumber >= 2){
    for(i in 2:colNumber){
      colCombinations = combinations(colnames(conditionsColumns[,1:colNumber]), i)
      for(j in 1:nrow(colCombinations)){
        interactionName = paste(colCombinations[j,], collapse = "_")
        designTable[,interactionName] = unite(conditionsColumns[,colCombinations[j,]], "merged", sep = ", ")$merged
      }
    }
    return(designTable)
  }
  else{
    return(designTable)
  }
}




# Method to sort count data and infoData according to the experimental setup and gff-file
sortThatData = function(rawCounts, infoData, gffData, gffType){
  # Purpose of this function is to sort the info data 
  # and set the column names of the raw data so they 
  # match with the info data (e.g. remove .bam ending)  
  
  # Additionally, row.names of count data will be set to gene name according to .gff-file
  # Additionally, sample prep table will receive a column with merged conditions
  
  #vec = numeric(0) #Character vector containing the correct order of column names. Will be used to sort info data 
  if(length(intersect(infoData[,1],colnames(rawCounts)))==0){
  for(i in 1:ncol(rawCounts)){                   # outer loop sets column of raw data... 
    checkCol = colnames(rawCounts)[i] 
    for(j in 1:nrow(infoData)){                  # ...which will be compared with nested loop using grepl
      if(grepl(infoData[j, 1], checkCol)){       # QBiC Code must be the first column of the info data for this to work
        #vec[i] = j                              # vec will contain the correct order the info data must be sorted with
        colnames(rawCounts)[i] = infoData[j, 1]  # replace column names with their correpsonding QBiC Code
      }
    }
  }
  }

  row.names(infoData) = infoData[,1]                  # set row names of info data to QBiC Code so it can be sorted by column names of count data
  # remove not required columns
  nonRequired = c("Chr", "Start", "End",	"Strand",	"Length",	"gene_name")
  rawCounts = rawCounts[,!colnames(rawCounts) %in% nonRequired]      
  print(rawCounts)
  print(infoData)
  #print(c("Geneid", row.names(infoData)))
  #infoData = infoData[colnames(rawCounts)[-1],]       # sort info data according to column name occurence in the counts file. Not occuring names will be removed
  infoData = infoData[row.names(infoData) %in% colnames(rawCounts)[-1],]       # remove all samples from sample prep file which are not in counts
  rawCounts = rawCounts[,c("Geneid", row.names(infoData))]
  
  # Change row names of raw counts to "locus_tag, gene name" (only tag, if gff-file has no corresponding gene): 
  if(gffType == "Bacteria"){
  names = gffData[gffData$locus_tag %in% rawCounts$Geneid & gffData$gbkey %in% c("Gene", "gene"),]$Name  # Match locus_tag of gff with Geneid and get gene names
  } else{
    names = gffData[gffData$gene_id %in% rawCounts$Geneid & gffData$gbkey %in% c("Gene", "gene"),]$gene_name  # Match locus_tag of gff with Geneid and get gene names
    
  }
  names[which(is.na(names))] = ""                                                                        # Otherwise the corresponding entry migh be NA (depends on gfffile) and might cause problems
  row.names(rawCounts) = rawCounts$Geneid
  if(length(names) > 0){
    # (if no matches between GeneID and gff$locus_tag were found, length(names) is 0) 
    combined_names = unite(data.frame(row.names(rawCounts), names), col = "merged", sep = ", ")  # merge tags and genenames (=> comma-seperated)
    row.names(rawCounts)[rawCounts$Geneid != names] = combined_names[rawCounts$Geneid != names,]
    row.names(rawCounts) = gsub("\\, $", "", rownames(rawCounts))  # remove commas if they are the last character
  }
  # Get condition-columns:
  is_treatment = grepl('condition', colnames(infoData), ignore.case = TRUE)
  treatments = data.frame(infoData[,is_treatment])
  
  # Change sample (column) names of rawCounts to sample prep (with suffix _1, _2, ... for replicates) and sort columns of counts according to info data
  merged_treatments = unite(treatments, "merged", sep = "_")
  merged_treatments = make.unique.2(merged_treatments$merged, sep = "_") # enumerates duplicates
  colnames(rawCounts)[-1] = merged_treatments
  
  # Add interaction-columns to design table
  infoData = addInteractionColumns(treatments, infoData)
  row.names(infoData) = merged_treatments

  return(list(rawCounts, infoData))     
}


# Method to split row index of a DF by comma into two columns => GeneID, Gene name
splitRowIndex = function(countTable){
  # get tags and gene names
  locusTags = gsub(',.*', '', row.names(rawCounts))
  geneNames = gsub('.*, ', '', row.names(rawCounts))
  # merge with countTable, remove potential additional geneid column w. tags
  newCounts = cbind(locusTags, geneNames, countTable[,!(tolower(colnames(countTable)) %in% "geneid")])
  return(newCounts)
}


# TPM normalization
normalizeTPM = function(rawCounts, gffFile, gffType){
  # normalize for gene length
  rpkTable = normalizeRPK(rawCounts,gffFile, gffType)
  
  # normalize for read depth
  totalSampleReadsPerMillion = colSums(rpkTable)/1e6 
  tpmTable = sweep(rpkTable, MARGIN = 2, totalSampleReadsPerMillion, FUN = "/")
  #tpmTable = tpmTable + 1  # add pseudo-counts
  return(tpmTable)
}

# Average coverage
normalizeRPK = function(rawCounts, gffFile, gffType){
  # remove locus tags that are not in rawCounts & calculate gene length
  if(gffType=="Bacteria"){
    gffFile = gffFile[!duplicated(gffFile$locus_tag),c("locus_tag", "start", "end")]
    idcol <- "locus_tag"
  }
  else{
    gffFile = gffFile[!duplicated(gffFile$gene_id),c("gene_id", "start", "end")]
    idcol <- "gene_id"
  }
  gffFile = gffFile[gffFile[,idcol] %in% rawCounts$Geneid,] 
  gffFile$length = gffFile$end - gffFile$start
  
  # normalize for gene length
  rpkTable = rawCounts[,-1]
  rpkTable = (rpkTable*1000)/gffFile$length
  return(rpkTable)
}

# Quantile normalization
# https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation = function(df){
  df_rank = apply(df,2,rank,ties.method="min")
  df_sorted = data.frame(apply(df, 2, sort))
  df_mean = apply(df_sorted, 1, mean)
  
  index_to_mean = function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final = apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) = rownames(df)
  return(df_final)
}


# Method to log-transform AND remove rows containing -Inf values 
logTransform = function(dataset){
  log2normCounts = log2(dataset)                   # produces -inf counts that were filtered (set to 0) by DESeq => set values to NA and ignore.
  log2normCounts[log2normCounts == -Inf] = NA      
  
  return(log2normCounts)
}


## RESULTS TABLES MANIPULATION ##

# Add new column "Gene name" to table (assumes that row.names is currently: 'locus tag, gene name'): 
addGeneNameCol = function(dds_results){
  if(dds_results[1,1] == "no significant genes found!"){
    return(dds_results)
  }
  else{
    # Separate row.names (usually they are composed by locus_tag, gene_name):
    splitVector <<- strsplit(row.names(dds_results), ", ")
    splitData = t(as.data.frame(splitVector))
    # Assign locus_tag ar rownames and add gene names:
    row.names(dds_results) = splitData[,1]
    if(ncol(splitData) >= 2){
      dds_results$'Gene name' = splitData[,2]           # sets gene name either to locus tag (if no gene has been found) or corresponding gene name in gff
    }
    else{
      dds_results$'Gene name' = row.names(dds_results)  # set gene name to locus tag, if no genes have been found
    }
    return(dds_results)
  }
}


# Add new column "Description" to table (by parsing 'product' of a given gff-file)
addDescriptionCol = function(dds_results, gff){
  if(dds_results[1,1] == "no significant genes found!"){
    return(dds_results)
  }
  else{
    # get CDS entries of gff-file and remove (potentially) duplicated locus tags
    gff = gff[!gff$gbkey %in% c("Gene", "gene"),]
    gff = gff[!duplicated(gff$locus_tag),]
    
    dds_results$Description = NA
    dds_results[row.names(dds_results) %in% gff$locus_tag,]$Description = gff[gff$locus_tag %in% row.names(dds_results),]$product
    return(dds_results)
  }
}


# Add new column "foldchange" to table
addFoldChangeCol = function(dds_results){
  dds_results$FoldChange = NA
  dds_results[dds_results$log2FoldChange >= 0 & !is.na(dds_results$log2FoldChange),]$FoldChange = 2^(dds_results[dds_results$log2FoldChange >= 0 & !is.na(dds_results$log2FoldChange),]$log2FoldChange)
  dds_results[dds_results$log2FoldChange < 0 & !is.na(dds_results$log2FoldChange),]$FoldChange = -2^abs(dds_results[dds_results$log2FoldChange < 0 & !is.na(dds_results$log2FoldChange),]$log2FoldChange)
  
  return(dds_results)
}


# Add new columsn for average TPM and log2 of average TPM. Only uses columns of specified contrast for calculation 
addAverageTPM = function(dds_results, tpmTable, contrast1, contrast2){
  # get matching columns, ignore spaces, commas and underscores to resolve potential ambiguities 
  tpmColumns = gsub(" |,|_", "", colnames(tpmTable))
  matchingCols = grepl(gsub(" |,|_", "", contrast1), tpmColumns) | grepl(gsub(" |,|_", "", contrast2), tpmColumns)
  
  dds_results$avgTPM = rowMeans(tpmTable[, matchingCols])
  return(dds_results)
}


# Function to sum up all results-table-manipulation functions and resort table:
extendAndSortResults = function(resultsData, gffFile, tpmData, contrast1, contrast2){
  resultsData = addGeneNameCol(resultsData)
  resultsData = addDescriptionCol(resultsData, gffFile)  
  resultsData = addFoldChangeCol(resultsData)
  resultsData = addAverageTPM(resultsData, tpmData, contrast1, contrast2)
  resultsData = resultsData[,c(7,1,2,9,3,10,4:6,8)]
  return(resultsData)
}


# Method to filter results data so it only contains significant genes (log FC >= 1 & p < alpha):
filterSignificantGenes = function(dds_results, alpha, logFCThreshold){
  
  dataset = as.data.frame(dds_results[!is.na(dds_results$padj),])
  significant_data = dataset[(abs(dataset$log2FoldChange) > logFCThreshold & dataset$padj < alpha), ]
  
  return(significant_data)
}


## OVERVIEW TABLE MANIPULATION ##

# Based on a (filtered) results-dataframe, make overview over up- and downregualted genes (single row of a dataframe):
significantOverview = function(dds_results, contrastVariable1, contrastVariable2){
  comparisonString = paste(contrastVariable1, "VS", contrastVariable2) 
  if(nrow(dds_results) == 0){
    up = 0
    down = 0
    total = 0
  }
  else{
    dataset = as.data.frame(dds_results)
    # Get infos from data:
    up =  nrow(dataset[(dataset$log2FoldChange > 0),])
    down = nrow(dataset[(dataset$log2FoldChange < 0),])
    total = nrow(dataset)
  }
  
  overview = data.frame(comparisonString, up, down, total)
  colnames(overview) = c("Conditions/Comparison", "UP", "DOWN", "TOTAL")
  return(overview)
}


# Update overview reactive value - objects: Append significantOverview object and update gene- and delete-button:
updateOverviewDataTable = function(overviewReactiveValues, significantOverviewEntry){
  overviewReactiveValues$data = rbind(overviewReactiveValues$data[,-c(1,6,7)], significantOverviewEntry)
  overviewReactiveValues$data$Set = make.unique.3(rep("Set", nrow(overviewReactiveValues$data)), sep = "_")
  overviewReactiveValues$data = overviewReactiveValues$data[,c(5,1:4)]
  overviewReactiveValues$data$Genes = shinyInput(actionButton, nrow(overviewReactiveValues$data), 'button_', label = "Show Genes", onclick = 'Shiny.setInputValue(\"genes_button\",  this.id.concat(\"_\", Math.random()))')
  overviewReactiveValues$data$Delete = shinyInput(actionButton, nrow(overviewReactiveValues$data), 'button_', label = "Delete", onclick = 'Shiny.setInputValue(\"delete_button\",  this.id.concat(\"_\", Math.random()))')
}


## PLOT RELATED ##

# Method that creates a volcano plot using ggplot2 
erupt = function(dds_results, logFCthreshold, alpha){
  # dds_results: results file from DESeq2
  # logFCthreshold: Specified by slide bar
  # alpha: Significance level, specified by slide bar
  # returns: Volcanoplot object and modified results table, that has everything that's required to make a volcano plot 
  
  fc_bound = logFCthreshold
  # transform results data
  res = as.data.frame(dds_results[!is.na(dds_results$padj),])
  # significance
  res$Expression = "NS"                                          
  is_up = res$log2FoldChange >= fc_bound & res$padj < alpha       
  is_down = res$log2FoldChange <= -fc_bound & res$padj < alpha    
  
  if(length(res[is_up,]$Expression)){
    res[is_up,]$Expression = "UP"          
  }
  if(length(res[is_down,]$Expression)){
    res[is_down,]$Expression = "DOWN"       
  }
  
  # data.frame for color palette. Will be checked with factor levels of res$expression 
  sign_colors = c("blue", "black", "red")
  sign_status = c("DOWN", "NS", "UP")
  col_assign = data.frame(sign_colors, sign_status)
  volc_palette = col_assign[col_assign$sign_status %in% levels(as.factor(res$Expression)),]$sign_colors
  
  # new column with -log10(padj) - interactive shiny CAN'T handle if you change parameters in aes() of ggplot and behaves weird
  res$neglog10_p_value = -log10(res$padj) 
  res$Locus_tag=row.names(res)
  res$Name=res[,"Gene name"]
  #print(res)
  # plot:
  vp = ggplot(res, aes(label= Locus_tag, label2=Name, x = log2FoldChange, y = neglog10_p_value, color = Expression, label3 = padj)) +
    # scatter:
    geom_point() +
    # lines (alpha and logFC):
    geom_hline(yintercept = (log10(alpha)/log10(10))*(-1), color = "darkgrey") +    # horizontal for alpha
    geom_vline(xintercept = c(fc_bound, -fc_bound), color = "darkgrey") +           # vertical for fc and -fc
    # color:
    scale_color_manual(values = volc_palette) +
    # axis formats:
    scale_x_continuous(breaks = c(round(min(res$log2FoldChange))):round(max(res$log2FoldChange))) # integers from rounded minimum to maximum of the log2FC
  return(list(vp, res))
}


# Create PCA plot (w. ggplot2) based on a pca dataset and a vector of variables with either one (color) or two variables (color+shape)
makePCA = function(pcaData, pcaGroups, fontSize = 11){
  # if the first variable is numeric, ggplot makes a color scale, which is not desired:
  pcaData[["data"]][[pcaGroups[1]]] = as.factor(pcaData[["data"]][[pcaGroups[1]]])
  
  # colors only, if there is only one group of interest:
  if(length(pcaGroups) == 1){
    pcaPlot = ggplot(pcaData[["data"]], aes(x = PC1, y = PC2, color = pcaData[["data"]][[pcaGroups[1]]])) +
      geom_point(size = 2) +
      labs(x = pcaData[["labels"]][["x"]], y = pcaData[["labels"]]["y"], color = pcaGroups) + 
      theme_grey(base_size = fontSize)
  }
  # colors + shapes, if there are two groups:
  else{
    # by default, ggplot runs out of shapes for >6 groups. Parsing to factor helps:
    pcaData[["data"]][[pcaGroups[2]]] = as.factor(pcaData[["data"]][[pcaGroups[2]]])
    
    pcaPlot = ggplot(pcaData[["data"]], aes(x = PC1, y = PC2, color = pcaData[["data"]][[pcaGroups[1]]], shape = pcaData[["data"]][[pcaGroups[2]]])) +
      geom_point(size = 2) +
      labs(x = pcaData[["labels"]][["x"]], y = pcaData[["labels"]]["y"], color = pcaGroups[1], shape = pcaGroups[2]) +
      scale_shape_manual(values = 1:nlevels(pcaData[["data"]][[pcaGroups[2]]])) +
      theme_grey(base_size = fontSize)
  }
  return(pcaPlot)
}


# Create venn diagram from list of gene-vectors and overview table (for labels)
makeVenn = function(geneList, overviewTable, fontsize = 11){
  vennDiagram = ggVennDiagram(geneList, 
                              #category.names = overviewTable$'Conditions/Comparison',  # CURRENTLY PROBLEMATIC WITH LONG NAMES
                              label = "count", set_size = fontsize, label_size = fontsize-5) + 
    theme(legend.position = "none")
  return(vennDiagram)
}


# Function to make people really upset
# Upset plot based on genelist and overview table (for labels)
makeUpset = function(geneList, overviewTable, fontsize = 11){
  names(geneList) = overviewTable$'Conditions/Comparison' #overviewTable[,1] 
  upsetPlot = upset(fromList(geneList), nsets = length(geneList), text.scale = fontsize/8)
  return(upsetPlot)
}


# Profile Plots (returns a colored plot and a black plot with mean line):
makeProfilePlots = function(tpmTable, geneList, summarize.replicates = TRUE, errorbars = FALSE, fontsize = 11){
  
  ## Data handling: 
  # Normalize
  tpmTable = log2(tpmTable)
  tpmTable = as.data.frame(quantile_normalisation(tpmTable))
  
  # Select specified genes and melt table
  tpmTable$Gene = gsub(".*, ", "", row.names(tpmTable))
  selectedTpm = tpmTable[tpmTable$Gene %in% geneList,]
  meanPerExperiment = as.data.frame(colMeans(selectedTpm[,1:ncol(selectedTpm)-1])) # for centroid (is basically already molten)
  selectedTpm = melt(selectedTpm, id.vars = "Gene")
  selectedTpm$Status = "__Gene"
  
  # Prepare centroid data set (must have same colnames and -order as selectedTPM)
  meanPerExperiment$Gene = "__Centroid"
  meanPerExperiment$Status = "Centroid"
  meanPerExperiment$variable = row.names(meanPerExperiment)
  meanPerExperiment = meanPerExperiment[,c(2, 4, 1, 3)]
  colnames(meanPerExperiment) = c("Gene", "variable", "value", "Status")
  
  if(summarize.replicates){
    # remove replicate suffix, which enables using the mean for gg-lineplots
    selectedTpm$variable = gsub("*_.$", "", selectedTpm$variable)
    meanPerExperiment$variable = gsub("*_.$", "", meanPerExperiment$variable)
  }
  
  ## Gene-wise (colored) plot:
  profilePlotColored = ggplot(data = selectedTpm, 
                       aes(x = variable, y = value, color = Gene, group = Gene)) +
    theme_grey(base_size = fontsize) +
    stat_summary(geom = 'line', fun = 'mean') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Quantile normalized log(TPM)") + xlab("Sample Condition") + ggtitle("Gene-wise expression profiles")
    
  
  ## Mean-expression plot:
  selectedTPMWithCentroid = rbind(selectedTpm, meanPerExperiment)
  profilePlotMean = ggplot(data = selectedTPMWithCentroid,
                              aes(x = variable, y = value, group = Gene, color = Status)) + 
    theme_grey(base_size = fontsize) +
    stat_summary(aes(size = Status, alpha = Status), geom = 'line', fun = 'mean') +
    scale_alpha_discrete(range = c(0.2, 1)) + 
    scale_size_discrete(range = c(0.4, 1.1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(values = c("black", "red")) + 
    ylab("Quantile normalized log(TPM)") + xlab("Sample Condition") + ggtitle("Mean of selected expression profiles") + 
    theme(legend.position = "none") 

  
  # Add errorbars (colored plot only)
  if(errorbars){
    profilePlotColored = profilePlotColored + stat_summary(geom = 'errorbar', fun.data = "mean_se", aes(width = 0.2))
  }
  
  # # Potential ToDo: Legend of for colored plot is cut off if > 14 genes are selected
  # if(length(geneList) > 14){
  #   profilePlotColored = profilePlotColored + theme(legend.position = "bottom")
  # }
  
  return(list(profilePlotColored, profilePlotMean))
}











