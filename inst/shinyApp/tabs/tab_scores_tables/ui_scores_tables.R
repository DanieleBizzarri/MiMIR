tab_scores_tables <- tabItem(
  tabName = "scores_tables",
  align = "center",
  
  #New analysis instructions
  h1("Metabolic Scores"),
  HTML('<hr style="border-color: #0088cc;">'),
  # fluidRow(
  #   column(
  #     width = 6,
  #     tabsetPanel(
  #       tabPanel(
  #         "Metabolic mortality score",
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         DT::dataTableOutput("mortScore_table") %>% withSpinner(),
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         style = "position: center; "
  #         )
  #       )
  #     ),
  #   column(
  #     width = 6,
  #     tabsetPanel(
  #       tabPanel(
  #         "MetaboAge",
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         DT::dataTableOutput("MetaboAge_table") %>% withSpinner(),
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         style = "position: center"
  #       )
  #     )
  #   )
  #   # ,
  #   # column(
  #   #   width = 4,
  #   #   tabsetPanel(
  #   #     tabPanel(
  #   #       "Surrogates",
  #   #       HTML('<hr style="border-color: #0088cc;">'),
  #   #       DT::dataTableOutput("surrogates_table") %>% withSpinner(),
  #   #       HTML('<hr style="border-color: #0088cc;">')
  #   #     )
  #   #   )
  #   # ),
  # ),
  # div(id="surro_table",
  #     tabsetPanel(
  #       tabPanel(
  #         "Surrogates",
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         DT::dataTableOutput("surrogates_table") %>% withSpinner(),
  #         HTML('<hr style="border-color: #0088cc;">')
  #       )
  #     )
  # ),
  # div(id="calib_table",
  #     tabsetPanel(
  #       tabPanel(
  #         "Calibrated Surrogates",
  #         HTML('<hr style="border-color: #0088cc;">'),
  #         DT::dataTableOutput("calibrated_surro_table") %>% withSpinner(),
  #         HTML('<hr style="border-color: #0088cc;">')
  #       )
  #     )
  # )
  
  br(),
  div(id = "data tabs",
      tabsetPanel(
        tabPanel(
          "Metabolic mortality score",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("mortScore_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center;"
          ),
        tabPanel(
          "MetaboAge",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("MetaboAge_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
          ),
        tabPanel(
          "Surrogates",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("surrogates_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">')
          ),
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




