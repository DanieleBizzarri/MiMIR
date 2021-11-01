tab_accuracy <- tabItem(
  tabName = "accuracy",
  align = "center",
  
  tabsetPanel(
    #Tab panel with the ROC curves of the surrogates
    tabPanel(
      title = "Surrogates ROC curves",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        sidebarPanel(
          width = 3,
          #Selection of the surrogates
          checkboxGroupInput(
            inputId = "surroc",
            label="Surrogates ROC curves:",
            choices=c("Surrogate sex" = "sex",
                      "Surrogate diabetes" = "diabetes",
                      "Surrogate lipidmed" = "lipidmed",
                      "Surrogate blood_pressure_lowering_med" = "blood_pressure_lowering_med",
                      "Surrogate current_smoking" = "current_smoking",
                      "Surrogate metabolic_syndrome" = "metabolic_syndrome",
                      "Surrogate alcohol_consumption" = "alcohol_consumption",
                      "Surrogate high_age" = "high_age",
                      "Surrogate middle_age" = "middle_age",
                      "Surrogate low_age" = "low_age",
                      "Surrogate obesity" = "obesity",
                      "Surrogate high_hscrp" = "high_hscrp",
                      "Surrogate high_triglycerides" = "high_triglycerides",
                      "Surrogate high_ldl_chol" = "high_ldl_chol",
                      "Surrogate low_hdlchol" = "low_hdlchol",
                      "Surrogate high_totchol" = "high_totchol",
                      "Surrogate low_eGFR" = "low_eGFR",
                      "Surrogate low_wbc" = "low_wbc",
                      "Surrogate low_hgb" = "low_hgb"
            ),
            selected = "sex",
          ),
          style = "text-align: left;"
        ),
        mainPanel(
          width = 9,
          plotlyOutput("ROCs", height = "600px") %>% withSpinner()
        )
      ),
      helpText("ROC curves of the surrogates calculated with the available phenotypic variables.
               On the side you can find the AUC. You can chose which surrogates ROC to visualize."),
      HTML('<hr style="border-color: #0088cc;">')
    ),
    #Tab panel with the t-test made with the surrogates
    tabPanel(
      title = "Surrogates t-tests",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 10,
        plotlyOutput("ttest", height = "500px") %>% withSpinner()
      ),
      helpText("In this plot we split the surrogates based on their original values and
               we calculate a t-test to see if they are well splitted")
    ),
    #Tab panel comparing accuracies of the surrogates in the uploaded datast to the accuracies in the LOBOV
    tabPanel(
      title = "Surrogates LOBOV comparison",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("LOBOV_surro", height = "600px") %>% withSpinner()
      ),
      helpText(" Boxplots comparing accuracies of the surrogates in the uploaded datast to the accuracies in the LOBOV")
    )
  )
)