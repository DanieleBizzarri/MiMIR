tab_surro_tables <- tabItem(
  tabName = "surro_tables",
  align = "center",
  
  # Metabolic scores tables
  h1("Metabolic Surrogates"),
  HTML('<hr style="border-color: #0088cc;">'),
  br(),
  div(id = "data tabs",
      tabsetPanel(
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




