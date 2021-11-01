tab_accuracy_scores <- tabItem(
  tabName = "accuracy_scores",
  align = "center",
  
  tabsetPanel(
    #Tab panel comparing MetaboAge and the age of the individuals
    tabPanel(
      title = "MetaboAge/age Scatterplot",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("scatter_metaboage", height = "600px") %>% withSpinner()
      ),
      helpText("Scatterplot to see how accurate MetaboAge is on your dataset. On the x axis the real age, on the y axis the predicted age.")
    ),
    #Histogram of the mortality score divided for 3 age ranges
    tabPanel(
      title = "Mortality score histogram per age range",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("hist_mort", height = "600px") %>% withSpinner()
      ),
      helpText("The distribution of the mortality score divided in age ranges: age<45, 45>age<65 and age>65")
    ),
    # Boxplot dividing one predictor for a binary variable
    tabPanel(
      title = "Predictor divided by binary variable",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        # Selection of which predictor to visualize
        sidebarPanel(
          width = 3,
          bootstrapPage(
            uiOutput("pred_choices")
          ),
          bootstrapPage(
            uiOutput("pheno_choices")
          )
        ),
        mainPanel(
          width = 9,
          plotlyOutput("ttest_score", height = "500px") %>% withSpinner()
          #htmlOutput("NAMES"),
        )
      ),
      helpText("Boxplot dividing the predicted score based on a binary variable."),
      HTML('<hr style="border-color: #0088cc;">')
      ),
    # Kaplan Meier with an event
    tabPanel(
      title = "Kaplan Meier",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        # Selection of which predictor to visualize
        sidebarPanel(
          width = 3,
          bootstrapPage(
            uiOutput("pred_km")
          ),
          # Event name
          textInput("Eventname", label = h3("Name of the Event"), value = "Event"),
        ),
        mainPanel(
          width = 9,
          plotlyOutput("km_score", height = "500px") %>% withSpinner()
        )
      ),
      helpText("Kaplan Meier of a score of your choice divided in tertiles associated to an event. You can chose the score to plot and the Event name.
               This figure is available only if in your phenotypes file you have a column named 'Event'(with 1 if the event happened, 0 if it didn't happen before censoring) 
               and one column 'EventAge' (with the age of the subject at censoring)"),
      HTML('<hr style="border-color: #0088cc;">')
    ),
    # MetaboWAS
    tabPanel(
      title = "MetaboWAS",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        # Selection of which predictor to visualize
        sidebarPanel(
          width = 3,
          bootstrapPage(
            uiOutput("metWASvar")
          ),
          bootstrapPage(
            uiOutput("metWAScovar")
          ),
          style = "text-align: left;"
        ),
        mainPanel(
          width = 9,
          plotlyOutput("manhplot", height = "500px") %>% withSpinner(),
          "MetaboWAS table",
          DT::dataTableOutput("res_metWAS_table") %>% withSpinner(),
        )
      ),
      helpText("Manhattan plot"),
      HTML('<hr style="border-color: #0088cc;">')
    )
    )
  )