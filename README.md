
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MiMIR

[![R-CMD-check](https://github.com/DanieleBizzarri/MiMIR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DanieleBizzarri/MiMIR/actions/workflows/R-CMD-check.yaml)

MiMIR (Metabolomics-based Models for Imputing Risk), is a a unique
graphical user interface that provides an intuitive framework for ad-hoc
statistical analysis of 1H-NMR metabolomics by Nightingale Health. It
allows to easily explore new metabolomics measurements assayed by
Nightingale Health; project previously published metabolic scores; and
calibrate the metabolic surrogate values to a desired dataset.

<img src="./inst/shinyApp/www/scaled_mimir_logo.svg" width="250" height="250" align="right">

To have a detail description of all the possible analyses available in
MiMIR, please take a look at the
Manual:<https://github.com/DanieleBizzarri/MiMIR/blob/main/man/MANUAL.pdf>
Please refer to our manuscripts when using these metabolic biomarkers in
your works: - mortality score: J. Deelen et al., ‘A metabolic profile of
all-cause mortality risk identified in an observational study of 44,168
individuals’, Nat. Commun., vol. 10, no. 1, pp. 1–8, Aug. 2019, doi:
10.1038/s41467-019-11311-9 - MetaboAge: van den Akker Erik B. et al.,
‘Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as
Biomarker of Age-related Disease’, Circ. Genomic Precis. Med., vol. 13,
no. 5, pp. 541–547, Oct. 2020, doi: 10.1161/CIRCGEN.119.002610. -
surrogate clinical variables: D. Bizzarri, M. J. T. Reinders, M.
Beekman, P. E. Slagboom, Bbmri-nl, and E. B. van den Akker, ‘1H-NMR
metabolomics-based surrogates to impute common clinical risk factors and
endpoints’, EBioMedicine, vol. 75, p. 103764, Jan. 2022, doi:
10.1016/j.ebiom.2021.103764. - COVID-severity score: Nightingale Health
UK Biobank Initiative, H. Julkunen, A. Cichońska, P. E. Slagboom, and P.
Würtz, ‘Metabolic biomarker profiling for identification of
susceptibility to severe pneumonia and COVID-19 in the general
population’, eLife, vol. 10, p. e63033, May 2021, doi:
10.7554/eLife.63033. - Type-2 diabetes score: A. V. Ahola-Olli et al.,
‘Circulating metabolites and the risk of type 2 diabetes: a prospective
study of 11,896 young adults from four Finnish cohorts’, Diabetologia,
vol. 62, no. 12, pp. 2298–2309, 2019, doi: 10.1007/s00125-019-05001-w. -
Cardiovascular event risk score: P. Würtz et al., ‘Metabolite profiling
and cardiovascular event risk: a prospective study of 3 population-based
cohorts’, Circulation, vol. 131, no. 9, pp. 774–785, Mar. 2015, doi:
10.1161/CIRCULATIONAHA.114.013116.

## Intalling

1.  Install the “devtools” package (if not already done):  

<!-- -->

    install.packages("devtools")

2.  Install the “MetaboRiSc” package:

<!-- -->

    library("devtools")
    devtools::install_github("DanieleBizzarri/MiMIR")

3.  Launch the application:

<!-- -->

    library("MetaboRiSc")
    MetaboRiSc::launchApp()

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

If you have problems in installing the applicationn, you can try
installing these packages manually:

    ## Shiny environment
    if (!require("shiny")) install.packages("shiny")
    if (!require("shinydashboard")) install.packages("shinydashboard")
    if (!require("shinyWidgets")) install.packages("shinyWidgets")
    if (!require("shinycssloaders")) install.packages("shinycssloaders")
    if (!require("shinyjs")) install.packages("shinyjs")

    ## Statistics libraries
    if (!require("DT")) install.packages("DT")
    if (!require("foreach")) install.packages("foreach")
    if (!require("glmnet")) install.packages("glmnet")
    if (!require("matrixStats")) install.packages("matrixStats")
    if (!require("plyr")) install.packages("plyr")
    if (!require("stats")) install.packages("stats")
    if (!require("reshape2")) install.packages("reshape2")
    if (!require("caret")) install.packages("caret")
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
