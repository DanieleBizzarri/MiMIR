# MiMIR

![image](./inst/shinyApp/www/MiMIR_logo.svg)

This package contains an Rshiny webtool developed to allow the
calculation of the metabolic predictorsdeveloped by the groups of MOLEPI
and LCBC (LUMC), from raw Nightingale Health 1H-NMR metabolomics data.
Please refer to our manuscripts when using these metabolic biomarkers in
your works: - mortality score: J. Deelen et al., ‘A metabolic profile of
all-cause mortality risk identified in an observational study of 44,168
individuals’, Nat. Commun., vol. 10, no. 1, pp. 1–8, Aug. 2019, doi:
10.1038/s41467-019-11311-9 - MetaboAge: van den Akker Erik B. et al.,
‘Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as
Biomarker of Age-related Disease’, Circ. Genomic Precis. Med., vol. 13,
no. 5, pp. 541–547, Oct. 2020, doi: 10.1161/CIRCGEN.119.002610. -
surrogate clinical variables: unpublished

## Intalling

1.  Install the “devtools” package (if not already done):  

<!-- end list -->

    install.packages("devtools")

2.  Install the “MiMIR” package:

<!-- end list -->

    library("devtools")
    devtools::install_github("DanieleBizzarri/MiMIR")

3.  Launch the application:

<!-- end list -->

    library("MiMIR")
    MiMIR::startApp()

## Quick Start

Note: By pressing the button “Dowload example” you can download a .zip
file, containing 2 files: the metabolic synthetic dataset, the
phenotypic synthetic dataset. These example dataset can be used to test
the App and to understand how the variables in your own dataset should
be named.

1.  Start the application
2.  Upload your metabolites with the same column names as in the example
    dataset (both CSV and TSV are accepted).
3.  Check if the App could find all the necessary metabolites in your
    dataset.
4.  Check if your dataset was correctly uploaded
5.  View the Predicted Scores and the Figures
6.  Download the results

## Requirements

R version: 3.6+

## Install packages

If you have problems in installing the application, you can try
installing these packages manually:

    ## Shiny environment
    if (!require("shiny")) install.packages("shiny")
    if (!require("shinydashboard")) install.packages("shinydashboard")
    if (!require("shinyWidgets")) install.packages("shinyWidgets")
    if (!require("shinycssloaders")) install.packages("shinycssloaders")
    if (!require("shinyjs")) install.packages("shinyjs")

    #Statistics libraries
    if (!require("DT")) install.packages("DT")
    if (!require("foreach")) install.packages("foreach")
    if (!require("glmnet")) install.packages("glmnet")
    if (!require("matrixStats")) install.packages("matrixStats")
    if (!require("plyr")) install.packages("plyr")
    if (!require("stats")) install.packages("stats")
    if (!require("reshape2")) install.packages("reshape2")
    if (!require("caret")) install.packages("caret")
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!require("limma")) BiocManager::install("limma", force=TRUE)
    if (!require("purrr")) install.packages("purrr")
    if (!require("dplyr")) install.packages("dplyr")
    if (!require("rmarkdown")) install.packages("rmarkdown")

    #Imaging libraries
    if (!require("pheatmap")) install.packages("pheatmap")
    if (!require("RColorBrewer")) install.packages("RColorBrewer")
    if (!require("pROC")) install.packages("pROC")
    if (!require("plotly")) install.packages("plotly")
    if (!require("heatmaply")) install.packages("heatmaply")
    if (!require("ggplot2")) install.packages("ggplot2")
    if (!require("ggfortify")) install.packages("ggfortify")
    if (!require("survival")) install.packages("survival")
    if (!require("survminer")) install.packages("survminer")
