#### =========== SHINYSEQ FUNCTIONS ========== ####
# Functions to be used in the server component

## MISC ##

# Altered version of make.unique
make.unique.2 = function(x, sep='.'){
  # https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
  # Makes unique names for duplicate entries in a vector. Default make.unique leaves the first duplicate unchanged. This function starts enumarting from duplicate 1. 
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
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

# Method to sort count data and infoData according to the experimental setup and gff-file
sortThatData = function(rawCounts, infoData, gffData){
  # Purpose of this function is to sort the info data 
  # and set the column names of the raw data so they 
  # match with the info data (e.g. remove .bam ending)  
  
  # Additionally, row.names of count data will be set to gene name according to .gff-file
  # Additionally, sample prep table will receive a column with merged conditions
  
  #vec = numeric(0) #Character vector containing the correct order of column names. Will be used to sort info data 
  
  for(i in 1:ncol(rawCounts)){                   # outer loop sets column of raw data... 
    checkCol = colnames(rawCounts)[i] 
    for(j in 1:nrow(infoData)){                  # ...which will be compared with nested loop using grepl
      if(grepl(infoData[j, 1], checkCol)){       # QBiC Code must be the first column of the info data for this to work
        #vec[i] = j                              # vec will contain the correct order the info data must be sorted with
        colnames(rawCounts)[i] = infoData[j, 1]  # replace column names with their correpsonding QBiC Code
      }
    }
  }
  
  row.names(infoData) = infoData[,1]                  # set row names of info data to QBiC Code so it can be sorted by column names of count data
  # row.names(rawCounts) = rawCounts[, 1]              # set row names of raw data to gene ID
  rawCounts = rawCounts[,-(2:7),]                     
  
  #print(c("Geneid", row.names(infoData)))
  #infoData = infoData[colnames(rawCounts)[-1],]       # sort info data according to column name occurence in the counts file. Not occuring names will be removed
  infoData = infoData[row.names(infoData) %in% colnames(rawCounts)[-1],]       # remove all samples from sample prep file which are not in counts
  rawCounts = rawCounts[,c("Geneid", row.names(infoData))]
  
  # Change row names of raw counts to "locus_tag, gene name" (only tag, if gff-file has no corresponding gene): 
  names = gffData[gffData$locus_tag %in% rawCounts$Geneid & gffData$gbkey == "Gene",]$Name  # Match locus_tag of gff with Geneid and get gene names
  row.names(rawCounts) = rawCounts$Geneid
  if(length(names) > 0){
    # (if no matches between GeneID and gff$locus_tag were found, length(names) is 0) 
    combined_names = unite(data.frame(row.names(rawCounts), names), col = "merged", sep = ", ")
    row.names(rawCounts)[rawCounts$Geneid != names] = combined_names[rawCounts$Geneid != names,]
  }
  
  
  # Add a merged-treatment column to sample prep:
  is_treatment = grepl('condition', colnames(infoData), ignore.case = TRUE)
  treatments = data.frame(infoData[,is_treatment])
  merged_treatments = unite(treatments, "merged", sep = ", ")
  infoData$All_conditions = merged_treatments$merged
  
  # Change sample (column) names of rawCounts to sample prep (with suffix _1, _2, ... for replicates) and sort columns of counts according to info data
  merged_treatments = unite(treatments, "merged", sep = "_")
  merged_treatments = make.unique.2(merged_treatments$merged, sep = "_") # enumarate duplicates
  colnames(rawCounts)[-1] = merged_treatments
  #rawCounts = rawCounts[,c("Geneid", merged_treatments)]  # sort 
  row.names(infoData) = merged_treatments
  
  return(list(rawCounts, infoData))     
}


# TPM normalization
normalizeTPM = function(rawCounts, gffFile){
  # remove locus tags that are not in rawCounts & calculate gene length
  gffFile = gffFile[!duplicated(gffFile$locus_tag),c("locus_tag", "start", "end")]
  gffFile = gffFile[gffFile$locus_tag %in% rawCounts$Geneid,] 
  gffFile$length = gffFile$end - gffFile$start
  
  # normalize for gene length
  rpkTable = rawCounts[,-1]
  rpkTable = (rpkTable*1000)/gffFile$length
  
  # normalize for read depth
  totalSampleReadsPerMillion = colSums(rpkTable)/1e6 
  tpmTable = sweep(rpkTable, MARGIN = 2, totalSampleReadsPerMillion, FUN = "/")
  #tpmTable = tpmTable + 1  # add pseudo-counts
  
  return(tpmTable)
}


# Method to log-transform AND remove rows containing -Inf values 
logTransform = function(dataset){
  log2normCounts = log2(dataset)                   # produces -inf counts that were filtered (set to 0) by DESeq => set values to NA and ignore.
  log2normCounts[log2normCounts == -Inf] = NA      
  
  return(log2normCounts)
}


## RESULTS TABLES MANIPULATION ##

# Method to filter results data so it only contains significant genes (log FC >= 1 & p < alpha):
filterSignificantGenes = function(dds_results, alpha, logFCThreshold){
  dataset = na.omit(dds_results)
  significant_data = dataset[(abs(dataset$log2FoldChange) > logFCThreshold & dataset$padj < alpha), ]
  
  if(length(significant_data)){
    return(significant_data)
  }
  else{
    return(data.frame("no significant genes found!"))
  }
}


# Based on a (filtered) results-dataframe, make overview over up- and downregualted genes (single row of a dataframe):
significantOverview = function(dds_results, contrastVariable1, contrastVariable2){
  comparisonString = paste(contrastVariable1, "VS", contrastVariable2) 
  
  if(dds_results[1,1] == "no significant genes found!"){
    up = 0
    down = 0
    total = 0
  }
  else{
    dataset = as.data.frame(na.omit(dds_results))
    # Get infos from data:
    up =  nrow(dataset[(dataset$log2FoldChange > 0),])
    down = nrow(dataset[(dataset$log2FoldChange < 0),])
    total = nrow(dataset)
  }
  
  overview = data.frame(comparisonString, up, down, total)
  colnames(overview) = c("Conditions/Comparison", "UP", "DOWN", "TOTAL")
  return(overview)
}


# Add new column "Gene name" to table (assumes that row.names is currently: 'locus tag, gene name'): 
addGeneNameCol = function(dds_results){
  if(dds_results[1,1] == "no significant genes found!"){
    return(dds_results)
  }
  else{
    # Separate row.names:
    splitVector = strsplit(row.names(dds_results), ", ")
    splitData = t(as.data.frame(splitVector))
    # Assign:
    row.names(dds_results) = splitData[,1]
    dds_results$'Gene name' = splitData[,2]
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
    gff = gff[gff$gbkey != "Gene",]
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

# Add new columsn for average TPM and log2 of average TPM 
addAverageTPM = function(dds_results, tpmTable){
  dds_results$avgTPM = rowMeans(tpmTable)
  #dds_results$'log2(avgTPM)' = log2(dds_results$avgTPM)
  return(dds_results)
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
  res = na.omit(as.data.frame(dds_results))
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
  sign_colors = c("red", "black", "blue")
  sign_status = c("DOWN", "NS", "UP")
  col_assign = data.frame(sign_colors, sign_status)
  volc_palette = col_assign[col_assign$sign_status %in% levels(as.factor(res$Expression)),]$sign_colors
  
  # new column with -log10(padj) - interactive shiny CAN'T handle if you change parameters in aes() of ggplot and behaves weird
  res$neglog10_p_value = -log10(res$padj) 
  # plot:
  vp = ggplot(res, aes(x = log2FoldChange, y = neglog10_p_value, color = Expression, tooltip = padj)) +
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
makePCA = function(pcaData, pcaGroups){
  # if the first variable is numeric, ggplot makes a color scale, which is not desired:
  pcaData[["data"]][[pcaGroups[1]]] = as.factor(pcaData[["data"]][[pcaGroups[1]]])
  
  # colors only, if there is only one group of interest:
  if(length(pcaGroups) == 1){
    pcaPlot = ggplot(pcaData[["data"]], aes(x = PC1, y = PC2, color = pcaData[["data"]][[pcaGroups[1]]])) +
      geom_point(size = 2) +
      theme(legend.title = element_blank()) + 
      labs(x = pcaData[["labels"]][["x"]], y = pcaData[["labels"]]["y"])
  }
  # colors + shapes, if there are two groups:
  else{
    pcaPlot = ggplot(pcaData[["data"]], aes(x = PC1, y = PC2, color = pcaData[["data"]][[pcaGroups[1]]], shape = pcaData[["data"]][[pcaGroups[2]]])) +
      geom_point(size = 2) +
      theme(legend.title = element_blank()) + 
      labs(x = pcaData[["labels"]][["x"]], y = pcaData[["labels"]]["y"])
  }
  return(pcaPlot)
}


# Create venn diagram from list of gene-vectors and overview table (for labels)
makeVenn = function(geneList, overviewTable){
  vennDiagram = ggVennDiagram(geneList, 
                              #category.names = overviewTable[,1],  # CURRENTLY PROBLEMATIC WITH LONG NAMES
                              label = "count") + 
    theme(legend.position = "none") 
  
  return(vennDiagram)
}

# Function to make people really upset
# Upset plot based on genelist and overview table (for labels)
makeUpset = function(geneList, overviewTable){
  names(geneList) = overviewTable[,1] 
  upsetPlot = upset(fromList(geneList), nsets = length(geneList))
  return(upsetPlot)
}











