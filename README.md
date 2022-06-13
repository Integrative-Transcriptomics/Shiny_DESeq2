# DESeq2-Vis
User-friendly shiny application for the interactive utilization of <a href=https://doi.org/10.1186/s13059-014-0550-8>DESeq2</a>. In addition to providing the standard DESeq2 methods, DESeq2-Vis provides further normalization options, as well as the visualization of gene profiles.

# Tutorial
The following section provides a quick rundown of the standard worklow and functions of DESeq2-Vis. When running DESeq2-Vis for the first time, all required packages will be installed automatically.

## 1. Data Upload & Analysis Paramenters

### 1.1 Upload your data
The following files are required for running DESeq2-Vis on your experimental data:
<ol>
  <li> A <tt>.tsv</tt>-file containing the raw counts for each sample in one column and locus tags as rownames. Tables acquired from <tt>featureCounts</tt> can also be directly uploaded to DESeq2-Vis.
  <li> An experimental design table (<tt>.tsv</tt>-file) containing sample names as row names and experimental conditions in columns. Each column containing an experimental condition to be analyzed must be indicated by containing the keyword <i>condition</i>. DESeq2-Vis will automatically scan this file and replace the column name of each sample in the counts-table by a merged string of the corresponding experimental conditions. A suffix indicating the replicate number is added to each sample name. Rows corresponding to samples that are not contained in the counts-table will be removed from the design table.
  <li> A GFF-file containing annotations and gene descriptions. Based on the locus tags, the corresponding gene name is added to each locus tag (unannotated genes will only contain the locus tag).
  <li> Press the <tt>Upload!</tt> button to start the upload- and scanning process. 
</ol>

### 1.2 Run DESeq on your data
<ol>
  <li> Select an <tt>Experimental Variable</tt> of your choice as a basis for the normalization and differential expression analysis. 
  <li> Select a normalization method. 
  <li> Adjust the significance level.
  <li> Press the <tt> Run DESeq!</tt> button in order to apply normalization.
</ol>

## 2. Normalization Results and Data Analyis
The <tt>Normalization</tt>-tab provides an overview of the normalized data and general data analysis tools for quality control and gene profile analysis.

<ul>
  <li> <tt>Normalized Counts</tt>: Contains normalized counts based on the specified method as well as a TPM-normalized table.
  <li> <tt>Boxplots</tt>: Provides simple boxplots of the normalized counts for each sample for quality control purposes (WIP).
  <li> <tt>PCA</tt>: Dot shape and color can be adjusted based on up to two experimental conditions. Use sliders to adjust plot- and font-size. Plot the PCA by pressing the <tt>Refresh Plot!</tt>-button.
  <li> <tt>Heatmaps</tt>: Contains two different, UPGMA-clustered heatmaps:
    <ul>
      <li> Pairwise distance between samples (euclidean distance of log2-normalized counts).
      <li> Heatmap of genes with the highest variance. Use the slider to adjust the amount of displayed genes and plot size.
    </ul>
   <li> <tt>Profile Plots</tt>: Use to display the the gene expression of one or multiple genes per condition or individual sample (use the checkbox <tt>Average replicates</tt> to switch). Gene profiles can be displayed as gene-wise expression profile (tab <tt>Sample Profiles</tt>) or <tt>Mean Expression Profile</tt> over all selected genes.  
</ul>

## 3. Differential Expression Analysis
The calculation of differential expression is carried out under <tt>Differential Expression</tt>-tab. <br>

The differential expression between two experimental groups is calculated as follows:
<ol>
  <li> Select a first and second experimental group. The log2-foldchange is calculated as $log$(1st group) - $log$(2nd group). 
  <li> Upon pressing <tt>Add to table!</tt>, the amount of signficantly up- and down-regulated genes as well as the total amount of significantly differentially expressed genes is display in tabular format. 
  <li> Press <tt>Show Genes</tt> to view the differential expression between the two conditions, gene descriptions as well as an interactive volcano plot.  
</ol>
Two or more differential expression results can be compared using a <tt>Venn Diagram</tt> or <tt>UpSet Plot</tt>. Sliders can be used to adjust the figure and font sizes.       


