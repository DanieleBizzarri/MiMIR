tab_distributions <- tabItem(
  tabName = "distributions",
  align = "center",
  
  tabsetPanel(
    # Tab Panel for the correlation of the predicted values
    tabPanel(
      title = "Predictor's correlations",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotlyOutput("heat_predictors", height = "600px") %>% withSpinner()
      )
    ),
    # Tab Panel for the missingness in the predicted values
    tabPanel(
      title = "Predictors' missingness",
      HTML('<hr style="border-color: #0088cc;">'),
      mainPanel(
        width = 9,
        plotOutput("heat_na_pred", height = "600px") %>% withSpinner()
      ),
      helpText("The missingness in the predicted values"),
      HTML('<hr style="border-color: #0088cc;">')
    ),
    # Tab Panel for the histograms of the predicted values
    tabPanel(
      title = "Predictors' histograms",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        # Selection of which predictor to visualize
        sidebarPanel(
          width = 3,
          checkboxGroupInput(
            inputId = "predictors",
            label="Metabolic predictors' histograms:",
            choices=c("Mortality score" = "mortScore",
              "MetaboAge" = "MetaboAge",
              "Surrogate sex" = "s_sex",
              "Surrogate diabetes" = "s_diabetes",
              "Surrogate lipidmed" = "s_lipidmed",
              "Surrogate blood_pressure_lowering_med" = "s_blood_pressure_lowering_med",
              "Surrogate current_smoking" = "s_current_smoking",
              "Surrogate metabolic_syndrome" = "s_metabolic_syndrome",
              "Surrogate alcohol_consumption" = "s_alcohol_consumption",
              "Surrogate high_age" = "s_high_age",
              "Surrogate middle_age" = "s_middle_age",
              "Surrogate low_age" = "s_low_age",
              "Surrogate obesity" = "s_obesity",
              "Surrogate high_hscrp" = "s_high_hscrp",
              "Surrogate high_triglycerides" = "s_high_triglycerides",
              "Surrogate high_ldl_chol" = "s_high_ldl_chol",
              "Surrogate low_hdlchol" = "s_low_hdlchol",
              "Surrogate high_totchol" = "s_high_totchol",
              "Surrogate low_eGFR" = "s_low_eGFR",
              "Surrogate low_wbc" = "s_low_wbc",
              "Surrogate low_hgb" = "s_low_hgb"
              ),
            selected = "mortScore",
          ),
          # Radio button to deciide if scale or not the scores
          radioButtons(
            inputId = "scaling",
            label="Do you want to see the scores scaled?",
            choices=c("Not scaled" = FALSE,
                      "Scaled" = TRUE
            ),
            selected = FALSE,
          ),
          
          style = "text-align: left;"
        ),
        mainPanel(
          width = 9,
          plotlyOutput("hist_predictors", height = "500px") %>% withSpinner()
        )
      ),
      helpText("Histograms to show the distributions of the predicted values. 
               You can look at each variables separately or together. 
               If you will look at them together you might want to look at them scaled, so they will have a similar range."),
      HTML('<hr style="border-color: #0088cc;">'),
      )
    )
)