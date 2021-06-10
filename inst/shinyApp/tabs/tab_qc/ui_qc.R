tab_qc <- tabItem(
  tabName = "qc",
  align = "center",
  
  #Pre-processing instructions
  h1("Pre-Processing"),
  uiOutput("qc_intro"),
  HTML('<hr style="border-color: #0088cc;">'),
  
  fluidRow(
    column(
      width = 6,
      h4("Quality Control for MetaboAge and the surrogates"),
      #QC of MetaboAge
      #Chose the threshold of zeros allowed
      sliderInput("Nmax_zero_metaboAge",
        label = "Set the number of missing values allowed per sample (default=Nmiss=1)",
        value = 1,
        min = 0,
        max = 55,
        step = 1
      ),
      #Chose the threshold of missing values allowed
      sliderInput("Nmax_miss_metaboAge",
        label = "Set the number of zeros allowed per sample (default=Nmiss=1)",
        value = 1,
        min = 0,
        max = 55,
        step = 1
      ),
      div(br(),
          htmlOutput("MetaboAge_settings"),
          br(),
          style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
      ),
      br(),
      #Print the selected thresholds and resulting dataset
      verbatimTextOutput("QC_metaboAge_text"),
      #style = "text-align:center"
    ),
    column(
      width = 6,
      h4("Quality Control for the Surrogates"),
      #QC of the surrogates
      #Chose the threshold of zeros allowed
      sliderInput("Nmax_zero_surrogates",
                  label = "Set the number of missing values allowed per sample (default=Nmiss=1)",
                  value = 1,
                  min = 1,
                  max = 55,
                  step = 1
      ),
      #Chose the threshold of missing values allowed
      sliderInput("Nmax_miss_surrogates",
                  label = "Set the number of zeros allowed per sample (default=Nmiss=1)",
                  value = 1,
                  min = 1,
                  max = 55,
                  step = 1
      ),
      div(br(),
          htmlOutput("Surrogates_settings"),
          br(),
          style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
      ),
      br(),
      #Print the selected thresholds and resulting dataset
      verbatimTextOutput("QC_surrogates_text"),
      #style = "text-align:center"
    )
  ),
  br(),
  HTML('<hr style="border-color: #0088cc;">'),
  actionButton(inputId = "run_button",
               label = "Run Analysis!"),
  br(),
  br()
)