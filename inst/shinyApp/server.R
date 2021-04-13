options(warn = -1)

server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
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
    scores_tables = "Metabolic scores",
    distributions = "Distributions",
    accuracy = "Accuracy",
    calibration = "Calibration",
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
          "Binary Phenotypes",
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
        "Predicted Scores",
        icon = icon("folder"),
        menuItem(
            "Metabolic scores",
            tabName = "scores_tables",
            icon = icon("table")
            ),
        menuItem(
          "Distributions",
          tabName = "distributions",
          icon = icon("chart-bar")
        ),
        menuItem(
          "Accuracy",
          tabName = "accuracy",
          icon = icon("chart-area")
        ),
        menuItem(
          "Surrogates Calibration",
          tabName = "calibration",
          icon = icon("chart-bar")
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
               menuSubItem(text = "MetaboAge", href = "https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.119.002610"),
               menuSubItem(text = "Mortality score", href = "https://www.nature.com/articles/s41467-019-11311-9"),
               menuSubItem(text = "Surrogate scores", href = "NaN")
      )
    )
  })
    
  source("tabs/tab_upload/server_upload.R", local = TRUE)
  source("tabs/tab_metabolites/server_metabolites.R", local = TRUE)
  source("tabs/tab_phenotypes/server_phenotypes.R", local = TRUE)
  source("tabs/tab_qc/server_qc.R", local = TRUE)
  source("tabs/tab_scores_tables/server_scores_tables.R", local = TRUE)
  source("tabs/tab_distributions/server_distributions.R", local = TRUE)
  source("tabs/tab_accuracy/server_accuracy.R", local = TRUE)
  source("tabs/tab_calibration/server_calibration.R", local = TRUE)
  source("tabs/tab_download/server_download.R", local = TRUE)
  source("tabs/tab_about/server_about.R", local = TRUE)
  }