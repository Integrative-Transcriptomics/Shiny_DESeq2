#### ========== PACKAGES FOR SHINYSEQ =========== #####

# install packages only if not already installed
packages = c("ggplot2", "pheatmap", "shiny", "shinythemes","shinyBS", "BiocManager", "tidyr", "ggVennDiagram", "UpSetR", "DT", "arrangements")           # required packages
req_packages = packages[!(packages %in% installed.packages()[,"Package"])]                                               # req_packages contains list of packages that are not already installed
if(length(req_packages)){install.packages(req_packages)}                                                                 # install req_packages
# Same for Bioconductor packages
biopackages = c("DESeq2", "rtracklayer", "preprocessCore")
req_biopackages = biopackages[!(biopackages %in% installed.packages()[,"Package"])]
if(length(req_biopackages)){BiocManager::install(req_biopackages)}

# load
library(arrangements)
library(DESeq2)
library(DT)
library(ggplot2)
library(ggVennDiagram)
library(pheatmap)
library(preprocessCore)
library(rtracklayer)
library(shiny)
library(shinyBS)
library(shinythemes)
library(tidyr)
library(UpSetR)

options(repos = BiocManager::repositories()) # required for publishing the tool  