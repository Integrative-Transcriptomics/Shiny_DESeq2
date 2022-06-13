# DESeq2-Vis
User-friendly shiny application for the interactive utilization of <a href=https://doi.org/10.1186/s13059-014-0550-8>DESeq2</a>. In addition to providing the standard DESeq2 methods, DESeq2-Vis provides further normalization options, as well as the visualization of gene profiles.

# Tutorial
The following section provides a quick rundown of the standard worklow and functions of DESeq2-Vis. When running DESeq2-Vis for the first time, all required packages will be installed automatically.

## 1. Data Upload
The following files are required for running DESeq2-Vis on your experimental data:
<ol>
  <li> A <tt>.tsv</tt>-file containing the raw counts for each sample in one column and locus tags as rownames. Tables acquired from <tt>featureCounts</tt> can also be directly uploaded to DESeq2-Vis.
  <li> An experimental design table (<tt>.tsv</tt>-file) containing sample names as row names and experimental conditions in columns. Each column containing an experimental condition to be analyzed must be indicated by containing the keyword <i>condition</i>. DESeq2-Vis will automatically scan this file and replace the column name of each sample in the counts-table by a merged string of the corresponding experimental conditions. A suffic indicating the replicate number is added to each sample name. Rows corresponding to samples that are contained in the counts-table will be removed from the design table.
  <li> A GFF-file containing annotations and gene descriptions. Based on the locus tags, the corresponding gene name is added to each locus tag (unannotated genes will only contain the locus tag).
  <li> Press the <tt>Upload!</tt> button to start the upload- and scanning process. 
</ol>

## 2. Running DESeq on your data
Select an <tt>Experimental Variable</tt> of your choice as a basis for the normalization 

