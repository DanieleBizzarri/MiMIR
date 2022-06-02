tab_calibration <- tabItem(
  tabName = "calibration",
  align = "center",
  
  tabsetPanel(
    # Tab Panel for the calibration reliability and histogram
    tabPanel(
      title = "Metabolic surrogates calibration",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        #It allows the selection of the surrogates
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "surro_cal",
            label="Metabolic surrogates calibration:",
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
            )
            ),
          numericInput("Nbins",
                      label = "Select the bins to for the probabilities",
                      value = 10
          ),
          style = "text-align: left;"
        ),
        mainPanel(
          width = 9,
          fluidPage(
            fluidRow(
              #Reliability plot of the calibrations
              loading_spin(plotlyOutput("reliability_calib"))
              ),
            fluidRow(
              #Histogram of the calibrations
              loading_spin(plotlyOutput("hist_calib"))
              #addSpinner(plotlyOutput("hist_calib"), spin = "circle", color = "#E41A1C")
              )
            )
        )
      ),
      helpText("Platt Calibration of the Surrogate values. 
               On the top you can see the Reliability diagram, on the bottom the histogram of the calibrations. \
               You can change the number of bins with the slider."),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
    # Tab Panel for the correlation of the calibrated surrogates
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
    # Tab Panel for the missingness of the calibrated surrogates
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