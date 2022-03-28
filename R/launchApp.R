## ----- START THE APPLICATION -----

#' Start MiMIR application.
#' 
#' @param launch.browser TRUE/FALSE
#' @return Opens application. If \code{launch.browser}=TRUE in default web browser
#'
#' @export

startApp <- function(launch.browser = TRUE) {
  appDir <- system.file("shinyApp", package = "MiMIR")
  
  if (appDir == "") {
    stop("Could not find 'MiMIR'. Try re-installing 'MiMIR'.",
         call. = FALSE)
  }
  
  library("MiMIR")
  options(shiny.maxRequestSize=30*1024^2)
  options(warn = -1)
  
  message("Initializing MiMIR...")
  suppressMessages(source(system.file("shinyApp/lib/libraries.R", package = "MiMIR")))
  
  if(launch.browser){
    shiny::runApp(
      appDir = appDir,
      host = "0.0.0.0",
      port = 1402,
      launch.browser = launch.browser
    )
  }else{
    shiny::runApp(
      appDir = appDir,
      host = "0.0.0.0",
      port = 1402,
      launch.browser = launch.browser
    )
  }
}