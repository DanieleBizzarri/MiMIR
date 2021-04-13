observe({
  checkReload()
  if (exists("predictors")) {
    shinyjs::enable("downloadTSV")
    shinyjs::enable("downloadCSV")
  } else {
    shinyjs::disable("downloadTSV")
    shinyjs::disable("downloadCSV")
  }
})

## create filename and save data as CSV
output$downloadCSV <- downloadHandler(
  filename = function() {
    paste0(input$downloadname, ".csv")
  },
  content = function(file) {
    write.table(
      predictors(),
      file,
      row.names = TRUE,
      col.names = NA,
      sep = ","
    )
  }
)

## create filename and save data as TSV
output$downloadTSV <- downloadHandler(
  filename = function() {
    paste0(input$downloadname, ".tsv")
  },
  content = function(file) {
    write.table(
      predictors(),
      file,
      row.names = TRUE,
      col.names = NA,
      sep = "\t"
    )
  }
)


output[["download_html"]] <- downloadHandler(
  # For PDF output, change this to "report.pdf"
  filename = function() {
    paste0("Metabolic_predictors_", gsub("-", "", Sys.Date()),".html")
  },
  if(dim(phenotypes()[1])!=0){
  content = function(file) {
    shiny::withProgress(
      message = paste0("Preparing Analysis report"),
      value = 0,
      {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          tempReport <- file.path("markdown/Rshiny_markdown.Rmd")
          file.copy("Rshiny_markdown.Rmd", tempReport, overwrite = TRUE)
          
          # Set up parameters to pass to Rmd document
          params<-list(metabo_measures=metabo_measures(),
                       phenotypes=phenotypes(),
                       bin_phenotypes=bin_phenotypes(),
                       bin_pheno_available=bin_pheno_available(),
                       mort_score=mort_score(),
                       MetaboAge=MetaboAge(),
                       surrogates=surrogates(),
                       predictors=predictors(),
                       calibrations=calibrations(),
                       Nbins=input$Nbins,
                       Nmax_miss_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_zero_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_miss_surrogates=input$Nmax_miss_metaboAge,
                       Nmax_zero_surrogates=input$Nmax_miss_metaboAge
          )
          shiny::incProgress(5/10)
          Sys.sleep(1)
          rmarkdown::render(tempReport,
                            output_file = file,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
          shiny::incProgress(10/10)
      })
  }
}
)

output[["download_html_no_pheno"]] <- downloadHandler(
  # For PDF output, change this to "report.pdf"
  filename = function() {
    paste0("Metabolic_predictors_", gsub("-", "", Sys.Date()),".html")
  },
    content = function(file) {
      shiny::withProgress(
        message = paste0("Preparing Analysis report"),
        value = 0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          tempReport <- file.path("markdown/Rshiny_markdown_no_pheno.Rmd")
          file.copy("Rshiny_markdown_no_pheno.Rmd", tempReport, overwrite = TRUE)
          # Set up parameters to pass to Rmd document
          params<-list(metabo_measures=metabo_measures(),
                       mort_score=mort_score(),
                       MetaboAge=MetaboAge(),
                       surrogates=surrogates(),
                       predictors=predictors(),
                       Nmax_miss_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_zero_metaboAge=input$Nmax_miss_metaboAge,
                       Nmax_miss_surrogates=input$Nmax_miss_metaboAge,
                       Nmax_zero_surrogates=input$Nmax_miss_metaboAge
          )
          shiny::incProgress(5/10)
          Sys.sleep(1)
          rmarkdown::render(tempReport,
                            output_file = file,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
          shiny::incProgress(10/10)
        }
      )
    }
)
# }else{
#   output[["download_html"]] <- downloadHandler(
#     # For PDF output, change this to "report.pdf"
#     filename = function() {
#       paste0("Metabolic_predictors_", gsub("-", "", Sys.Date()),".html")
#     },
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
#                        )
#           shiny::incProgress(5/10)
#           Sys.sleep(1)
#           rmarkdown::render(tempReport,
#                     output_file = file,
#                     params = params,
#                     envir = new.env(parent = globalenv())
#                     )
#           shiny::incProgress(10/10)
#         }
#       )
#     })
# }

# output[["download_html"]] <- downloadHandler(
#   # For PDF output, change this to "report.pdf"
#   filename = "report.html",
#   content = function(file) {
#     # Copy the report file to a temporary directory before processing it, in
#     # case we don't have write permissions to the current working dir (which
#     # can happen when deployed).
#     tempReport <- file.path("markdown/report.Rmd")
#     file.copy("report.Rmd", tempReport, overwrite = TRUE)
# 
#     # Set up parameters to pass to Rmd document
#     params <- list(metabo_measures=metabo_measures())
# 
#     # Knit the document, passing in the `params` list, and eval it in a
#     # child of the global environment (this isolates the code in the document
#     # from the code in this app).
#     rmarkdown::render(tempReport, output_file = file,
#                       params = params,
#                       envir = new.env(parent = globalenv())
#     )
#   }
#   )

