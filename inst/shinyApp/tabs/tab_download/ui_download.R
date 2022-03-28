tab_download <- tabItem(
  tabName = "download",
  align = "center",
  
  HTML('<hr style="border-color: #0088cc;">'),
  h1("Download predictors tables"),
  shinySaveButton("download_pred", "Save predictors table", "Save file as...", 
                  filetype = list(csv = ".csv"), filename = paste0("metabolic_predictors_", Sys.Date()),
                  viewtype = "icon"),
  br(),
  br(),
  shinySaveButton("download_calib", "Save calibrated surrogates", "Save file as...", 
                  filetype = list(csv = ".csv"), filename = paste0("metabolic_calibrated_surrogates_", Sys.Date()),
                  viewtype = "icon"),
  HTML('<hr style="border-color: #0088cc;">'),
  h2("Download Report of the analyses computed"),
  br(),
  #Download the Rmarkdown if the user loaded the phenotypic file
  shinySaveButton("download_html", "Download Analysis Report with accuracies", "Save file as...", 
                  filetype = list(html = ".html"), filename = paste0("MiMIR_metabolic_predictors_report_", Sys.Date()),
                  viewtype = "icon"),
  br(),
  br(),
  #Download the Rmarkdown if the user didn't load the phenotypic file
  shinySaveButton("download_html_no_pheno", "Download Analysis Report with accuracies without accuracy analyses", "Save file as...", 
                  filetype = list(html = ".html"), filename = paste0("MiMIR_metabolic_predictors_report_", Sys.Date()),
                  viewtype = "icon"),
  br(),
  br(),
  HTML('<hr style="border-color: #0088cc;">')
)
