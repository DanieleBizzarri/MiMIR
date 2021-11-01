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
        # T2D AholaOlli tab
        tabPanel(
          "T2Dscore",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("T2Dscore_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
        )
        ,
        # CVD score tab
        tabPanel(
          "CVDscore",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("CVDscore_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
        ),
        # COVID score tab
        tabPanel(
          "COVIDscore",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("COVIDscore_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
        ),
        # All predictors tab
        tabPanel(
          "All predictors",
          HTML('<hr style="border-color: #0088cc;">'),
          DT::dataTableOutput("predictors_table") %>% withSpinner(),
          HTML('<hr style="border-color: #0088cc;">'),
          style = "position: center"
        )
        )
  )
)




