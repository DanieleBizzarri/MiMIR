# Enable or disable the predictor table button
observe({
  checkReload()
  if (required()) {
    shinyjs::enable("download_pred")
  } else {
    shinyjs::disable("download_pred")
    
  }
})

# Write the prediction table
shiny::observeEvent(
  eventExpr = input$download_pred,
  handlerExpr = {
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    session = getDefaultReactiveDomain()
    shinyFiles::shinyFileSave(input = input,id = "download_pred",roots = volumes,session = session,
                              restrictions = system.file(package = "base"))
      
    fileinfo <- shinyFiles::parseSavePath(roots = volumes, selection = input$download_pred)
      
      if (nrow(fileinfo) > 0){
        #if(input$downloadCSV2$filetype == ".csv"){
          write.table(predictors(), fileinfo$datapath, row.names = TRUE, col.names = NA, sep = ",")
        # }else{
        #   write.table(
        #     predictors(),
        #     file,
        #     row.names = TRUE,
        #     col.names = NA,
        #     sep = "\t"
        #   )
        # }
      }
  }
)

# Enable or disable the calibration table button
observe({
  checkReload()
  if (required() & !is.null(phenotypes())) {
    shinyjs::enable("download_calib")
  } else {
    shinyjs::disable("download_calib")
  }
})
# Write the calibration table
shiny::observeEvent(
  eventExpr = input$download_calib,
  handlerExpr = {
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    session = getDefaultReactiveDomain()
    shinyFiles::shinyFileSave(input = input,id = "download_calib",roots = volumes,session = session,
                              restrictions = system.file(package = "base"))
    fileinfo <- shinyFiles::parseSavePath(roots = volumes, selection = input$download_calib)
    if (nrow(fileinfo) > 0){
      write.table(calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available()),
        fileinfo$datapath, row.names = TRUE, col.names = NA,sep = ",")
    }
  }
)

# Enable or disable the Report button
observe({
  checkReload()
  if (required() & !is.null(phenotypes())) {
    shinyjs::enable("download_html")
  } else {
    shinyjs::disable("download_html")
  }
})
# Write the Report with phenotypes
shiny::observeEvent(
  eventExpr = input$download_html,
  handlerExpr = {
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    session = getDefaultReactiveDomain()
    shinyFiles::shinyFileSave(input = input,id = "download_html",roots = volumes,session = session,
                              restrictions = system.file(package = "base"))
    
    
    fileinfo <- shinyFiles::parseSavePath(roots = volumes, selection = input$download_html)
    if (nrow(fileinfo) > 0){
      shiny::withProgress(message = paste0("Preparing Analysis report"), value = 0, {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          tempReport <- file.path("markdown/Rshiny_markdown.Rmd")
          file.copy("Rshiny_markdown.Rmd", tempReport, overwrite = TRUE)
          # Set up parameters to pass to Rmd document
          params<-list(metabo_measures=metabo_measures(), phenotypes=phenotypes(), bin_phenotypes=bin_phenotypes(),
                       bin_pheno_available=bin_pheno_available(), mort_score=mort_score(), MetaboAge=MetaboAge(),
                       surrogates=surrogates(), predictors=predictors(), calibrations=calibrations(), Nbins=input$Nbins,
                       Nmax_miss_metaboAge=input$Nmax_miss_metaboAge, Nmax_zero_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_miss_surrogates=input$Nmax_miss_metaboAge,Nmax_zero_surrogates=input$Nmax_miss_metaboAge
          )
          shiny::incProgress(5/10)
          Sys.sleep(1)
          rmarkdown::render(tempReport,output_file = fileinfo$datapath, params = params,envir = new.env(parent = globalenv()))
          shiny::incProgress(10/10)
        })
    }
  }
)

# Enable or disable the Report button without phenotypes information
observe({
  checkReload()
  if (required()) {
    shinyjs::enable("download_html_no_pheno")
  } else {
    shinyjs::disable("download_html_no_pheno")
  }
})
#write Report button without phenotypes information
shiny::observeEvent(
  eventExpr = input$download_html_no_pheno,
  handlerExpr = {
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    session = getDefaultReactiveDomain()
    shinyFiles::shinyFileSave(input = input,id = "download_html_no_pheno",roots = volumes,session = session,
                              restrictions = system.file(package = "base"))
    
    
    fileinfo <- shinyFiles::parseSavePath(roots = volumes, selection = input$download_html_no_pheno)
    if (nrow(fileinfo) > 0){
      shiny::withProgress(message = paste0("Preparing Analysis report"), value = 0,{
          shiny::incProgress(1/10)
          Sys.sleep(1)
          tempReport <- file.path("markdown/Rshiny_markdown_no_pheno.Rmd")
          file.copy("Rshiny_markdown_no_pheno.Rmd", tempReport, overwrite = TRUE)
          # Set up parameters to pass to Rmd document
          params<-list(metabo_measures=metabo_measures(),mort_score=mort_score(),MetaboAge=MetaboAge(),
                       surrogates=surrogates(),predictors=predictors(), Nmax_miss_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_zero_metaboAge=input$Nmax_miss_metaboAge, Nmax_miss_surrogates=input$Nmax_miss_metaboAge, Nmax_zero_surrogates=input$Nmax_miss_metaboAge)
          shiny::incProgress(5/10)
          Sys.sleep(1)
          rmarkdown::render(tempReport, output_file = fileinfo$datapath, params = params, envir = new.env(parent = globalenv()))
          shiny::incProgress(10/10)
        })
    }
  }
)

# output$downloadCSV2 <-  downloadHandler(
#     filename = function() {
#       paste0(input$downloadCSV2, ".csv")
#     },
#     content = function(file) {
#       write.table(
#         predictors(),
#         file,
#         row.names = TRUE,
#         col.names = NA,
#         sep = ","
#       )
#     }
#   )

# ## output all the predictors and save data as CSV
# output$downloadCSV <- downloadHandler(
#   filename = function() {
#     paste0(input$downloadname, ".csv")
#   },
#   content = function(file) {
#     write.table(
#       predictors(),
#       file,
#       row.names = TRUE,
#       col.names = NA,
#       sep = ","
#     )
#   }
# )
# 
# ## output all the predictors and save data as TSV
# output$downloadTSV <- downloadHandler(
#   filename = function() {
#     paste0(input$downloadname, ".tsv")
#   },
#   content = function(file) {
#     write.table(
#       predictors(),
#       file,
#       row.names = TRUE,
#       col.names = NA,
#       sep = "\t"
#     )
#   }
# )
# 
# # Check if the metabolites are available
# observe({
#   checkReload()
#   if (required()) {
#     shinyjs::enable("downloadTSV_calib")
#     shinyjs::enable("downloadCSV_calib")
#   } else {
#     shinyjs::disable("downloadTSV_calib")
#     shinyjs::disable("downloadCSV_calib")
#   }
# })
# 
# ## Output the calibrated surrogates and save data as CSV
# output$downloadCSV_calib <- downloadHandler(
#   filename = function() {
#     paste0(input$downloadname_calib, ".csv")
#   },
#   content = function(file) {
#     write.table(
#       calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available()),
#       file,
#       row.names = TRUE,
#       col.names = NA,
#       sep = ","
#     )
#   }
# )
# 
# ## create the calibrated surrogates and save data as TSV
# output$downloadTSV_calib <- downloadHandler(
#   filename = function() {
#     paste0(input$downloadname_calib, ".tsv")
#   },
#   content = function(file) {
#     write.table(
#       calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available()),
#       file,
#       row.names = TRUE,
#       col.names = NA,
#       sep = "\t"
#     )
#   }
# )

# observe({
#   checkReload()
#   if (required() & !is.null(phenotypes())) {
#     shinyjs::enable("download_html")
#   } else {
#     shinyjs::disable("download_html")
#   }
# })
# 
# # Outputs the histogram of the 
# output[["download_html"]] <- downloadHandler(
#   # For PDF output, change this to "report.pdf"
#   filename = function() {
#     paste0("Metabolic_predictors_", gsub("-", "", Sys.Date()),".html")
#   },
#   content = function(file) {
#     shiny::withProgress(
#       message = paste0("Preparing Analysis report"),
#       value = 0,
#       {
#           shiny::incProgress(1/10)
#           Sys.sleep(1)
#           tempReport <- file.path("markdown/Rshiny_markdown.Rmd")
#           file.copy("Rshiny_markdown.Rmd", tempReport, overwrite = TRUE)
#           
#           # Set up parameters to pass to Rmd document
#           params<-list(metabo_measures=metabo_measures(),
#                        phenotypes=phenotypes(),
#                        bin_phenotypes=bin_phenotypes(),
#                        bin_pheno_available=bin_pheno_available(),
#                        mort_score=mort_score(),
#                        MetaboAge=MetaboAge(),
#                        surrogates=surrogates(),
#                        predictors=predictors(),
#                        calibrations=calibrations(),
#                        Nbins=input$Nbins,
#                        Nmax_miss_metaboAge=input$Nmax_miss_metaboAge,
#                        Nmax_zero_metaboAge=input$Nmax_miss_metaboAge,
#                        Nmax_miss_surrogates=input$Nmax_miss_metaboAge,
#                        Nmax_zero_surrogates=input$Nmax_miss_metaboAge
#           )
#           shiny::incProgress(5/10)
#           Sys.sleep(1)
#           rmarkdown::render(tempReport,
#                             output_file = file,
#                             params = params,
#                             envir = new.env(parent = globalenv())
#           )
#           shiny::incProgress(10/10)
#       })
# }
# )
# 
# observe({
#   checkReload()
#   if (required()) {
#     shinyjs::enable("download_html_no_pheno")
#   } else {
#     shinyjs::disable("download_html_no_pheno")
#   }
# })
# 
# output[["download_html_no_pheno"]] <- downloadHandler(
#   # For PDF output, change this to "report.pdf"
#   filename = function() {
#     paste0("Metabolic_predictors_", gsub("-", "", Sys.Date()),".html")
#   },
#     content = function(file) {
#       shiny::withProgress(
#         message = paste0("Preparing Analysis report"),
#         value = 0,
#         {
#           shiny::incProgress(1/10)
#           Sys.sleep(1)
#           tempReport <- file.path("markdown/Rshiny_markdown_no_pheno.Rmd")
#           file.copy("Rshiny_markdown_no_pheno.Rmd", tempReport, overwrite = TRUE)
#           # Set up parameters to pass to Rmd document
#           params<-list(metabo_measures=metabo_measures(),
#                        mort_score=mort_score(),
#                        MetaboAge=MetaboAge(),
#                        surrogates=surrogates(),
#                        predictors=predictors(),
#                        Nmax_miss_metaboAge=input$Nmax_miss_metaboAge,
#                        Nmax_zero_metaboAge=input$Nmax_miss_metaboAge,
#                        Nmax_miss_surrogates=input$Nmax_miss_metaboAge,
#                        Nmax_zero_surrogates=input$Nmax_miss_metaboAge
#           )
#           shiny::incProgress(5/10)
#           Sys.sleep(1)
#           rmarkdown::render(tempReport,
#                             output_file = file,
#                             params = params,
#                             envir = new.env(parent = globalenv())
#           )
#           shiny::incProgress(10/10)
#         }
#       )
#     }
# )
