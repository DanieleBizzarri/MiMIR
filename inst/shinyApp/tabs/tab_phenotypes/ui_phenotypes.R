tab_phenotypes <- tabItem(
  tabName = "phenotypes",
  align = "center",
  
  tabsetPanel(
    # Tab Panel for the summary of the phenotypic data
    tabPanel(
      title = "Phenotypes summary",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        DT::dataTableOutput("summ_pheno") %>% withSpinner()
      ),
      helpText("The summary of the phenotypic table you uploaded"),
      #HTML('<hr style="border-color: #0088cc;">')
    ),
    # Tab Panel for the binary phenotypes
    tabPanel(
      title = "Binarized phenotype barplot",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("bin_phen_barplot", height = "600px") %>% withSpinner()
      ),
      helpText("Binarized phenotypes barplots"),
      #HTML('<hr style="border-color: #0088cc;">')
    ),
    # Tab Panel for the correlation of the phenotypes uploaded
    tabPanel(
      title = "Binary phenotypes correlations",
      HTML('<hr style="border-color: #0088cc;">'),
        mainPanel(
          width = 9,
          plotlyOutput("heat_pheno", height = "600px") %>% withSpinner()
        ),
      helpText("The correlations of the phenotypes used in our prediction models"),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
    # Tab Panel for the missingness of the phenotypes uploaded
  tabPanel(
    title = "Binary phenotypes missing values",
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
