tab_download <- tabItem(
  tabName = "download",
  align = "center",
  
  HTML('<hr style="border-color: #0088cc;">'),
  h1("Download files"),
  br(),
  br(),
  #input for the name of the predictors table
  textInput("downloadname", "Input a name for the predictors table", value = paste0("metabolic_predictors_", Sys.Date())),
  br(),
  #Download the table in CSV or TSV
  downloadButton("downloadCSV", label = "Download predictors as CSV"),
  downloadButton("downloadTSV", label = "Download predictors as TSV"),
  br(),
  #input for the name of the calibrated surrogates
  textInput("downloadname_calib", "Input a name for the calibrated surrogates table", value = paste0("metabolic_calibrated_surrogates_", Sys.Date())),
  br(),
  #Download the table in CSV or TSV
  downloadButton("downloadCSV_calib", label = "Download calibrated surrogates as CSV"),
  downloadButton("downloadTSV_calib", label = "Download calibrated surrogates as TSV"),
  br(),
  br(),
  HTML('<hr style="border-color: #0088cc;">'),
  h2("Download R Markdown analysis"),
  br(),
  #Download the Rmarkdown if the user loaded the phenotypic file
  downloadButton("download_html", label = "Download Analysis Report with accuracy analyses"),
  br(),
  br(),
  #Download the Rmarkdown if the user didn't load the phenotypic file
  downloadButton("download_html_no_pheno", label = "Download Analysis Report without accuracy analyses"),
  br(),
  br(),
  HTML('<hr style="border-color: #0088cc;">')
)
