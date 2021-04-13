tab_heatmaps <- tabItem(
  tabName = "heatmaps",
  align = "center",
  
  
  tabsetPanel(
    tabPanel(
      title = "Predictor's correlations",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("heat_predictors", height = "600px") %>% withSpinner()
      )
    )
    # tabPanel(
    #   title = "Metabolites Correlations",
    #   HTML('<hr style="border-color: #0088cc;">'),
    #   sidebarLayout(
    #     position = "right",
    #     sidebarPanel(
    #       width = 3,
    #         radioButtons(
    #           "MET56_14_cor",
    #           "Metabolites selections:",
    #           c("MetaboAge and surrogates selection" = "MET56",
    #             "Mortality selection" = "MET14")
    #         )
    #       ),
    #   mainPanel(
    #     width = 9,
    #     plotlyOutput("heat_met", height = "600px") %>% withSpinner()
    #     )
    #   ),
    #   helpText("The correlations of the metabolites used in our prediction models, you can chose the set of the mortality score or the one of MetaboAge/Surrogates"),
    #   HTML('<hr style="border-color: #0088cc;">'),
    # ),
    # tabPanel(
    #   title = "Metabolites NAs",
    #   HTML('<hr style="border-color: #0088cc;">'),
    #   sidebarLayout(
    #     position = "right",
    #     sidebarPanel(
    #       width = 3,
    #       radioButtons(
    #         "MET56_14_nas",
    #         "Metabolites selections:",
    #         c("MetaboAge and surrogates selection" = "MET56",
    #           "Mortality selection" = "MET14")
    #         )
    #       ),
    #     mainPanel(
    #       width = 9,
    #       plotOutput("heat_na_metabo", height = "600px") %>% withSpinner()
    #       )
    #     ),
    #   helpText("The missingness of the metabolites used in our prediction models, you can chose the set of the mortality score or the one of MetaboAge/Surrogates"),
    #   HTML('<hr style="border-color: #0088cc;">')
    #   ),
    # tabPanel(
    #   title = "Phenotypes NAs",
    #   HTML('<hr style="border-color: #0088cc;">'),
    #   mainPanel(
    #     width = 9,
    #     plotOutput("heat_na_pheno", height = "600px") %>% withSpinner()
    #   )
    #   )
  )
)