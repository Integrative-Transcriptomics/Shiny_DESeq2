#### =========== SHINYSEQ FUNCTIONS ========== ####
# Functions to be used in the server component



# Alter version of make.unique
make.unique.2 = function(x, sep='.'){
  # https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
  # Makes unique names for duplicate entries in a vector. Default make.unique leaves the first duplicate unchanged. This function starts enumarting from duplicate 1. 
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}



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
  # row.names(rawCounts) = rawCounts[, 1]             # set row names of raw data to gene ID
  rawCounts = rawCounts[,-(2:7)]                     
  infoData = infoData[colnames(rawCounts)[-1],]       # sort info data according to column name occurence in the counts file. Not occuring names will be removed
  
  # Change row names of raw counts to "locus_tag, gene name" (only tag, if gff-file has no corresponding gene): 
  names = gffData[gffData$locus_tag %in% rawCounts$Geneid & gffData$gbkey == "Gene",]$Name  # Match locus_tag of gff with Geneid and get gene names
  row.names(rawCounts) = rawCounts$Geneid
  if(length(names) > 0){
    # (if no matches between GeneID and gff$locus_tag were found, length(names) is 0) 
    combined_names = unite(data.frame(row.names(rawCounts), names), col = "merged", sep = ", ")
    row.names(rawCounts)[rawCounts$Geneid != names] = combined_names[rawCounts$Geneid != names,]
  }
  
  
  # Add a merged-treatment column to sample prep:
  is_treatment = grepl('Condition..treatment', colnames(infoData))
  treatments = infoData[,is_treatment]
  merged_treatments = unite(treatments, "merged", sep = ", ")
  infoData$All_conditions = merged_treatments$merged
  
  # Change sample (column) names of rawCounts to sample prep (with suffix _1, _2, ... for replicates)
  merged_treatments = unite(treatments, "merged", sep = "_")
  merged_treatments = make.unique.2(merged_treatments$merged, sep = "_") # enumarate duplicates
  colnames(rawCounts)[-1] = merged_treatments
  row.names(infoData) = merged_treatments
  
  return(list(rawCounts, infoData))     
}



# Method to log-transform AND remove rows containing -Inf values 
logTransform = function(dataset){
  log2normCounts = log2(dataset)                   # produces -inf counts that were filtered (set to 0) by DESeq => set values to NA and ignore.
  log2normCounts[log2normCounts == -Inf] = NA      
  
  return(log2normCounts)
}


# Method to filter results data so it only contains significant genes (log FC >= 1 & p < alpha)
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