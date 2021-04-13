## ----- START THE APPLICATION -----


#' Start metaboRiSc application.

startApp <- function(launch.browser = TRUE) {
  appDir <- system.file("shinyApp", package = "metaboRiSc")
  
  if (appDir == "") {
    stop("Could not find 'metaboRiSc'. Try re-installing 'metaboRiSc'.",
         call. = FALSE)
  }
  
  library("metaboRiSc")
  #options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  #options(spinner.color = "#0088cc")
  options(warn = -1)
  
  message("Initializing metaboRiSc...")
  suppressMessages(source(system.file("shinyApp/lib/libraries.R", package = "metaboRiSc")))
  
  # load("inst/shinyApp/data/PARAM__2018-06-18_02-16-17.457.RData")
  # usethis::use_data(PARAM)
  # PARAM_surrogates<-readRDS("inst/shinyApp/data/PARAM_surrogates_2021_03_26.RData")
  # usethis::use_data(PARAM_surrogates)
  
  
  shiny::runApp(
    appDir = appDir,
    host = "0.0.0.0",
    port = 1402,
    launch.browser = launch.browser
  )
}