tab_scores_tables <- tabItem(
  tabName = "scores_tables",
  align = "center",
  
  # Metabolic scores tables
  h1("Metabolic Scores"),
  HTML('<hr style="border-color: #0088cc;">'),
  br(),
  div(id = "data tabs",
      tabsetPanel(
        #Mortality score tab
        tabPanel(
          "Metabolic mortality score",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("mortScore_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center;"
          ),
        # MetaboAge tab
        tabPanel(
          "MetaboAge",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("MetaboAge_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
          ),
        # Surrogates tab
        tabPanel(
          "Surrogates",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("surrogates_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">')
          ),
        # Calibrated surrogates tab
        tabPanel(
          "Calibrated Surrogates",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("calibrated_surro_table") %>% withSpinner(),
          helpText("These are the Metabolic Surrogate markers that were calibrated on this dataset, using the Platt Calibration"),
          HTML('<hr style="border-color: #0088cc;">')
          )
        )
  )
)




