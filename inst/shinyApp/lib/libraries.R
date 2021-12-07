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

