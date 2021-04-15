## ----- START THE APPLICATION -----


#' Start metaboRiSc application.

startApp <- function(launch.browser = TRUE) {
  appDir <- system.file("shinyApp", package = "MetaboRiSc")
  
  if (appDir == "") {
    stop("Could not find 'MetaboRiSc'. Try re-installing 'MetaboRiSc'.",
         call. = FALSE)
  }
  
  library("MetaboRiSc")
  #options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  options(warn = -1)
  
  message("Initializing metaboRiSc...")
  suppressMessages(source(system.file("shinyApp/lib/libraries.R", package = "MetaboRiSc")))
  
  # load("inst/shinyApp/data/PARAM__2018-06-18_02-16-17.457.RData")
  # usethis::use_data(PARAM, MetaboRiSc, overwrite=T)
  # PARAM_surrogates<-readRDS("inst/shinyApp/data/PARAM_surrogates_2021_04_13.RData")
  # usethis::use_data(PARAM_surrogates)
  
  
  shiny::runApp(
    appDir = appDir,
    host = "0.0.0.0",
    port = 1402,
    launch.browser = launch.browser
  )
}