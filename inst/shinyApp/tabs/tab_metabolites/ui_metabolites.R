tab_metabolites <- tabItem(
  tabName = "metabolites",
  align = "center",
  
  tabsetPanel(
    # Tab Panel for the correlation of the metabolite concentrations uploaded
    tabPanel(
      title = "Metabolites Correlations",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        sidebarPanel(
          width = 3,
          # Metabolites selection of MetaboAge and the surrogates or the one of the Mortality score
          radioButtons(
            "MET56_14_cor",
            "Metabolites selections:",
            c("MetaboAge and surrogates selection" = "MET56",
              "Mortality selection" = "MET14")
          )
        ),
        mainPanel(
          width = 9,
          plotlyOutput("heat_met", height = "600px") %>% withSpinner()
        )
      ),
      helpText("The correlations of the metabolites used in our prediction models, you can chose the set of the mortality score or the one of MetaboAge/Surrogates"),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
    # Tab Panel for the missingness of the metabolite concentrations uploaded
    tabPanel(
      title = "Metabolites NAs",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        sidebarPanel(
          width = 3,
          # Metabolites selection of MetaboAge and the surrogates or the one of the Mortality score
          radioButtons(
            "MET56_14_nas",
            "Metabolites selections:",
            c("MetaboAge and surrogates selection" = "MET56",
              "Mortality selection" = "MET14")
          )
        ),
        mainPanel(
          width = 9,
          plotOutput("heat_na_metabo", height = "600px") %>% withSpinner()
        )
      ),
      helpText("The missingness of the metabolites used in our prediction models, you can chose the set of the mortality score or the one of MetaboAge/Surrogates"),
      HTML('<hr style="border-color: #0088cc;">')
    ),
    # Tab Panel for the histograms of the metabolite concentrations uploaded
  tabPanel(
    title = "Metabolites histograms",
    HTML('<hr style="border-color: #0088cc;">'),
    sidebarLayout(
      position = "right",
      sidebarPanel(
        width = 3,
        #Radio button if you want to scale them or not
        radioButtons(
          inputId = "scale_met",
          label="Do you want to see the scores scaled?",
          choices=c("Not scaled" = FALSE,
                    "Scaled" = TRUE
          ),
          selected = FALSE,
        ),
        #Selection of which metabolites to print
        checkboxGroupInput(
          inputId = "metabolites",
          label="Metabolites' histograms:",
          choices=c("Total lipids in medium VLDL" = "m_vldl_l",
                    "Total lipids in small VLDL" = "s_vldl_l",
                    "Total lipids in very small VLDL" = "xs_vldl_l",
                    "Total lipids in IDL" = "idl_l",
                    "Total cholesterol in IDL" = "idl_c",
                    "Total lipids in large LDL" = "l_ldl_l",
                    "Total lipids in medium LDL" = "m_ldl_l",
                    "Total lipids in medium HDL" = "m_hdl_l",
                    "Total lipids in small HDL" = "s_hdl_l",
                    "Total lipids in small LDL" = "s_ldl_l",
                    "Mean diameter for VLDL particles" = "vldl_d",
                    "Mean diameter for LDL particles" = "ldl_d",
                    "Mean diameter for HDL particles" = "hdl_d",
                    "Serum total cholesterol" = "serum_c",
                    "Total cholesterol in VLDL" = "vldl_c",
                    "Total cholesterol in LDL" = "ldl_c",
                    "Total cholesterol in HDL" = "hdl_c",
                    "Total cholesterol in HDL2" = "hdl2_c",
                    "Total cholesterol in HDL3" = "hdl3_c",
                    "Total phosphoglycerides" = "totpg",
                    "Phosphatidylcholine and other cholines" = "pc",
                    "Sphingomyelins" = "sm",
                    "Total cholines" = "totcho",
                    "Apolipoprotein A-I" = "apoa1",
                    "Apolipoprotein B" = "apob",
                    "Total fatty acids" = "totfa",
                    "Estimated degree of unsaturation" = "unsatdeg",
                    "docosahexaenoic acid" = "dha",
                    "linoleic acid" = "la",
                    "Omega-3 fatty acids" = "faw3",
                    "Omega-6 fatty acids" = "faw6",
                    "Polyunsaturated fatty acids" = "ufa",
                    "Monounsaturated fatty acids" = "mufa",
                    "Saturated fatty acids" = "sfa",
                    "Ratio of omega-3 fatty acids to total fatty acids" = "faw3fa",
                    "Apolipoprotein B" = "faw6_fa",
                    "Ratio of omega-6 fatty acids to total fatty acids" = "totfa",
                    "Ratio of polyunsaturated fatty acids to total fatty acids" = "pufa_fa",
                    "Ratio of monounsaturated fatty acids to total fatty acids" = "mufa_fa",
                    "Ratio of saturated fatty acids to total fatty acids" = "sfa_fa",
                    "Glucose" = "glc",
                    "Lactate" = "lac",
                    "Citrate" = "cit",
                    "Alanine" = "ala",
                    "Glutamine" = "gln",
                    "Histidine" = "his",
                    "Isoleucine" = "ile",
                    "Leucine" = "leu",
                    "Valine" = "val",
                    "Phenylalanine" = "phe",
                    "Tyrosine" = "tyr",
                    "Acetate" = "ace",
                    "Acetoacetate" = "acace",
                    "Creatinine" = "crea",
                    "Albumin" = "alb",
                    "gp" = "Glycoprotein acetyls",
                    "Total lipids in chylomicrons and extremely large VLDL" = "xxl_vldl_l"
                    ),
          selected = "m_vldl_l",
          ),
        style = "text-align: left;"
      ),
      mainPanel(
        width = 9,
        plotlyOutput("hist_metabolites", height = "500px") %>% withSpinner()
      )
    ),
    helpText("Histograms to show the distributions of the metabolites uploaded. 
               You can also look at each variables together. 
               If you will look at them together you might want to look at them scaled, so they will have a similar range."),
    HTML('<hr style="border-color: #0088cc;">'),
  )
)
)