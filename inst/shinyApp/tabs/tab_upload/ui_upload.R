tab_upload <- tabItem(
  tabName = "upload",
  align = "center",
  
  # h1("Welcome"),
  # uiOutput("welcome_text"),
  
  #Upload files
  h1("Upload your files"),
  
  fluidPage(
    br(),
    column(
      width = 3,
      div(
        id = "upload_metabolites",
        br(),
        fileInput(
          inputId = "file_samples",
          label = "Choose your metabolic files:",
          multiple = FALSE,
          accept = c(
            "CSV",
            ".csv",
            "TSV",
            ".tsv",
            "TXT",
            ".txt",
            "ZIP",
            ".zip",
            "GZ",
            ".gz"
          )
        )
      ),
      div(
        id = "upload_phenotypes",
        br(),
        fileInput(
          inputId = "file_phenotypes",
          label = "Choose your phenotype file (optional):",
          multiple = FALSE,
          accept = c(
            "CSV",
            ".csv",
            "TSV",
            ".tsv",
            "TXT",
            ".txt",
            "ZIP",
            ".zip",
            "GZ",
            ".gz"
          )
        )
      ),
      br(),
      h5(strong("Example dataset:")),
      downloadButton('downloadData', 'Download example'),
      br(),
    ),
    column(
      width = 9,
      uiOutput("get_started"),
      br(),
    ),
    style = "position: center; border-radius: 20px; border: 2px solid #0088cc;'width:100px;'",
  ),
  HTML('<hr style="border-color: #0088cc;">'),
  fluidPage(
    h3(strong("Column names found")),
    textOutput("required_met"),
    br(),
    column(
      width = 6,
      div(
        DT::dataTableOutput("found_met") %>% withSpinner()
      )
      ),
    column(
      width = 6,
      div(
        DT::dataTableOutput("found_phen") %>% withSpinner(),
      )
    ),
    ),
  helpText("Red= are required to continue; Orange=  not found but optional; Green= found"),
  br(),
  div(id = "data tabs",
      tabsetPanel(
        tabPanel(
          "Metabolic measurements",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("metabo_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">')
        ),
        tabPanel(
          "Phenotypic variables",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("pheno_table") %>% withSpinner(),
          helpText("All the phenotypic variables available in your dataset"),
          HTML('<hr style="border-color: #0088cc;">')
        ),
        tabPanel(
          "Binarized phenotypic variables",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("bin_pheno_table") %>% withSpinner(),
          helpText("These is the set of the binarized phenotypic variables available in your dataset that we will use for the accuracy"),
          HTML('<hr style="border-color: #0088cc;">')
        )
        
      )
  )
  
  
  )
