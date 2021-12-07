options(warn = -1)

server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  #Check the reload of the files
  checkReload <- function() {
    
    is.null(input$file_samples)
    is.null(input$file_phenotypes)
    
  }
  
  ## All page names
  pages <- list(
    upload = "Data upload",
    metabolites= "Metabolites",
    phenotypes = "Phenotypes",
    qc = "Quality Control",
    surro_tables = "Metabolic surrogates",
    distributions = "Distributions",
    accuracy = "Accuracy",
    calibration = "Calibration",
    scores_tables = "Metabolic scores",
    distributions_scores = "Distributions",
    accuracy_scores = "Accuracy",
    download = "Download results",
    about = "About"
  )
  
  # Render current page name to ui
  output[["current_page"]] <- renderText({
    tryCatch({
      page_name <- pages[[input$sidebar]]
    }, error = function(err) {
      return(NULL)
    })
  })
  
  ## Render navigation bar to ui
  output[["sidebar_tabs"]] <- renderMenu({
    sidebarMenu(
      id = "sidebar",
      menuItem(
        "Data upload",
        tabName = "upload",
        icon = icon("upload")
      ),
      menuItem(
        "Your Dataset",
        icon = icon("folder"),
        menuItem(
          "Metabolites",
          tabName = "metabolites",
          icon = icon("chart-bar")
        ),
        menuItem(
          "Phenotypes",
          tabName = "phenotypes",
          icon = icon("chart-bar")
        )
      ),
      menuItem(
        "Settings QC (optional)",
        tabName = "qc",
        icon = icon("cogs")
      ),
      menuItem(
        "Metabolomics Surrogates",
        icon = icon("folder"),
        menuItem(
            "Metabolomics surrogates",
            tabName = "surro_tables",
            icon = icon("table")
            ),
        menuItem(
          "Distributions/Correlations",
          tabName = "distributions",
          icon = icon("chart-bar")
        ),
        menuItem(
          "Accuracy",
          tabName = "accuracy",
          icon = icon("chart-area")
        ),
        menuItem(
          "Calibration",
          tabName = "calibration",
          icon = icon("chart-bar")
        )
      ),
      menuItem(
        "Metabolomics Scores",
        icon = icon("folder"),
        menuItem(
          "Metabolomics scores",
          tabName = "scores_tables",
          icon = icon("table")
        ),
        menuItem(
          "Distributions/Correlations",
          tabName = "distributions_scores",
          icon = icon("chart-bar")
        ),
        menuItem(
          "Miscellaneus items",
          tabName = "accuracy_scores",
          icon = icon("chart-area")
        )
      ),
      menuItem(
        "Download Results",
        tabName = "download",
        icon = icon("download")
      ),
      menuItem(
        "About",
        tabName = "about",
        icon = icon("info-circle")
      ),
      menuItem("References",startExpanded = TRUE,
               tabName = "references",
               icon = icon("book-open"),
               menuSubItem(text = "Surrogate scores", href = "https://www.medrxiv.org/content/10.1101/2021.07.19.21258470v1"),
               menuSubItem(text = "MetaboAge", href = "https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.119.002610"),
               menuSubItem(text = "Mortality score", href = "https://www.nature.com/articles/s41467-019-11311-9"),
               menuSubItem(text = "T2D score", href = "https://link.springer.com/article/10.1007/s00125-019-05001-w"),
               menuSubItem(text = "CVD score", href = "https://www.ahajournals.org/doi/10.1161/circulationaha.114.013116"),
               menuSubItem(text = "COVID score", href = "https://elifesciences.org/articles/63033")
      )
    )
  })
    
  #Load all the server tabs
  source("tabs/tab_upload/server_upload.R", local = TRUE)
  source("tabs/tab_metabolites/server_metabolites.R", local = TRUE)
  source("tabs/tab_phenotypes/server_phenotypes.R", local = TRUE)
  source("tabs/tab_qc/server_qc.R", local = TRUE)
  source("tabs/tab_surro_tables/server_surro_tables.R", local = TRUE)
  source("tabs/tab_scores_tables/server_scores_tables.R", local = TRUE)
  source("tabs/tab_distributions/server_distributions.R", local = TRUE)
  source("tabs/tab_distributions_scores/server_distributions_scores.R", local = TRUE)
  source("tabs/tab_accuracy/server_accuracy.R", local = TRUE)
  source("tabs/tab_accuracy_scores/server_accuracy_scores.R", local = TRUE)
  source("tabs/tab_calibration/server_calibration.R", local = TRUE)
  source("tabs/tab_download/server_download.R", local = TRUE)
  source("tabs/tab_about/server_about.R", local = TRUE)
  }