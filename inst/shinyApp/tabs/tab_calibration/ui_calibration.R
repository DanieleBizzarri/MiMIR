tab_calibration <- tabItem(
  tabName = "calibration",
  align = "center",
  
  tabsetPanel(
    tabPanel(
      title = "Metabolic surrogates calibration",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        sidebarPanel(
          width = 3,
          radioButtons(
            inputId = "surrogates",
            label="Metabolic surrogates' calibration:",
            choices=c("Surrogate sex" = "s_sex",
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
            )
            ),
          sliderInput("Nbins",
                      label = "Select the bins to for the probabilities",
                      value = 10,
                      min = 5,
                      max = 100,
                      step = 1
          ),
          style = "text-align: left;"
        ),
        mainPanel(
          width = 9,
          #plotlyOutput("hist_predictors", height = "500px") %>% withSpinner()
          fluidPage(
            fluidRow(
              plotlyOutput("reliability_calib") %>% withSpinner()
              ),
            fluidRow(
              plotlyOutput("hist_calib") %>% withSpinner()
              )
            )
        )
      ),
      helpText("Platt Calibration of the Surrogate values. 
               On the top you can see the Reliability diagram, on the bottom the histogram of the calibrations. \
               You can change the number of bins with the slider."),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
      tabPanel(
        title = "Calibrated surrogates Correlations",
        HTML('<hr style="border-color: #0088cc;">'),
        mainPanel(
          width = 9,
          plotlyOutput("heat_calib", height = "600px") %>% withSpinner()
        ),
        helpText("The correlations of the calibrated surrogates"),
        HTML('<hr style="border-color: #0088cc;">'),
      ),
      tabPanel(
        title = "Calibrated surrogates NAs",
        HTML('<hr style="border-color: #0088cc;">'),
        mainPanel(
          width = 9,
          plotOutput("heat_na_calib", height = "600px") %>% withSpinner()
        ),
        helpText("The missingness in the calibrated surrogates"),
        HTML('<hr style="border-color: #0088cc;">')
      )
  )
  
)