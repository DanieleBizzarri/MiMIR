tab_phenotypes <- tabItem(
  tabName = "phenotypes",
  align = "center",
  
  tabsetPanel(
    tabPanel(
      title = "Phenotypes Correlations",
      HTML('<hr style="border-color: #0088cc;">'),
        mainPanel(
          width = 9,
          plotlyOutput("heat_pheno", height = "600px") %>% withSpinner()
        ),
      helpText("The correlations of the phenotypes used in our prediction models"),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
  tabPanel(
    title = "Phenotypes NAs",
    HTML('<hr style="border-color: #0088cc;">'),
    mainPanel(
      width = 9,
      plotOutput("heat_na_pheno", height = "600px") %>% withSpinner()
    ),
    helpText("The missingness in the phenotypic dataset you uploaded"),
    HTML('<hr style="border-color: #0088cc;">')
    )
    )
)