## ----- START THE APPLICATION -----

#' startMiMIR
#' 
#' Start the application MiMIR.
#' 
#' @param launch.browser TRUE/FALSE
#' @return Opens application. If \code{launch.browser}=TRUE in default web browser
#'
#' @details 
#' This function starts the R-Shiny tool called MiMIR (Metabolomics-based Models for Imputing Risk), 
#' a graphical user interface that provides an intuitive framework for ad-hoc statistical analysis of Nightingale Health's 1H-NMR metabolomics data
#' and allows for the projection and calibration of 24 pre-trained metabolomics-based models, without any pre-required programming knowledge. 
#'
#' @references 
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, doi: 10.1038/s41467-019-11311-9.
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, doi: 10.1007/s00125-019-05001-w
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, doi: 10.1161/CIRCULATIONAHA.114.013116
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, doi: 10.1016/j.ebiom.2021.103764
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, doi:10.1161/CIRCGEN.119.002610
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