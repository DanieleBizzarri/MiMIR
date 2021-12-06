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
            c("MetaboAge and surrogates" = "MET56",
              "Mortality score" = "MET14",
              "COVID score"="MET_COVID",
              "T2D score" = "MET_T2D",
              "CVD score" = "MET_CVD")
          )
        ),
        mainPanel(
          width = 9,
          plotlyOutput("heat_met", height = "600px") %>% withSpinner()
        )
      ),
      helpText("The correlations of the available metabolic features set for each metabolic score"),
      HTML('<hr style="border-color: #0088cc;">'),
    ),
    # Tab Panel for the missingness of the metabolite concentrations uploaded
    tabPanel(
      title = "Metabolites missing  values",
      HTML('<hr style="border-color: #0088cc;">'),
      sidebarLayout(
        position = "right",
        sidebarPanel(
          width = 3,
          # Metabolites selection of MetaboAge and the surrogates or the one of the Mortality score
          radioButtons(
            "MET56_14_nas",
            "Metabolites:",
            c("MetaboAge and surrogates" = "MET56",
              "Mortality" = "MET14",
              "COVID score "="MET_COVID",
              "T2D score" = "MET_T2D",
              "CVD score" = "MET_CVD")
          )
        ),
        mainPanel(
          width = 9,
          plotOutput("heat_na_metabo", height = "600px") %>% withSpinner()
        )
      ),
      helpText("The missingness in all the metabolites used for the scores"),
      HTML('<hr style="border-color: #0088cc;">')
    ),
    # Tab Panel for the histograms of the metabolite concentrations uploaded
  tabPanel(
    title = "Metabolites Distributions",
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
                    "scaled" = TRUE
          ),
          selected = FALSE,
        ),
        #Selection of which metabolites to print
        checkboxGroupInput(
          inputId = "metabolites",
          label="Metabolites' histograms:",
          choices=c("Total lipids in medium VLDL (m_vldl_l)" = "m_vldl_l",
                    "Total lipids in small VLDL (s_vldl_l)" = "s_vldl_l",
                    "Total lipids in very small VLDL (xs_vldl_l)" = "xs_vldl_l",
                    "Total lipids in IDL (idl_l)" = "idl_l",
                    "Total cholesterol in IDL (idl_c)" = "idl_c",
                    "Total lipids in large LDL (l_ldl_l)" = "l_ldl_l",
                    "Total lipids in medium LDL (m_ldl_l)" = "m_ldl_l",
                    "Total lipids in medium HDL (m_hdl_l)" = "m_hdl_l",
                    "Total lipids in small HDL (s_hdl_l)" = "s_hdl_l",
                    "Total lipids in small LDL (s_ldl_l)" = "s_ldl_l",
                    "Mean diameter for VLDL particles (vldl_d)" = "vldl_d",
                    "Mean diameter for LDL particles (ldl_d)" = "ldl_d",
                    "Mean diameter for HDL particles (hdl_d)" = "hdl_d",
                    "Serum total cholesterol (serum_c)" = "serum_c",
                    "Total cholesterol in VLDL (vldl_c)" = "vldl_c",
                    "Total cholesterol in LDL (ldl_c)" = "ldl_c",
                    "Total cholesterol in HDL (hdl_c)" = "hdl_c",
                    "Total cholesterol in HDL2 (hdl2_c)" = "hdl2_c",
                    "Total cholesterol in HDL3 (hdl3_c)" = "hdl3_c",
                    "Total phosphoglycerides (totpg)" = "totpg",
                    "Phosphatidylcholine and other cholines (pc)" = "pc",
                    "Sphingomyelins (sm)" = "sm",
                    "Total cholines (totcho)" = "totcho",
                    "Apolipoprotein A-I (apoa1)" = "apoa1",
                    "Apolipoprotein B (apob)" = "apob",
                    "Total fatty acids (totfa)" = "totfa",
                    "Estimated degree of unsaturation (unsatdeg)" = "unsatdeg",
                    "docosahexaenoic acid (dha)" = "dha",
                    "linoleic acid (la)" = "la",
                    "Omega-3 fatty acids (faw3)" = "faw3",
                    "Omega-6 fatty acids (faw6)" = "faw6",
                    "Polyunsaturated fatty acids (pufa)" = "pufa",
                    "Monounsaturated fatty acids (mufa)" = "mufa",
                    "Saturated fatty acids (sfa)" = "sfa",
                    "Ratio of omega-3 fatty acids to total fatty acids (faw3_fa)" = "faw3_fa",
                    "Ratio of omega-6 fatty acids to total fatty acids (totfa)" = "totfa",
                    "Ratio of polyunsaturated fatty acids to total fatty acids (pufa_fa)" = "pufa_fa",
                    "Ratio of monounsaturated fatty acids to total fatty acids (mufa_fa)" = "mufa_fa",
                    "Ratio of saturated fatty acids to total fatty acids (sfa_fa)" = "sfa_fa",
                    "Glucose (glc)" = "glc",
                    "Lactate (lac)" = "lac",
                    "Citrate (cit)" = "cit",
                    "Alanine (ala)" = "ala",
                    "Glutamine (gln)" = "gln",
                    "Histidine (his)" = "his",
                    "Isoleucine (ile)" = "ile",
                    "Leucine (leu)" = "leu",
                    "Valine (val)" = "val",
                    "Phenylalanine (phe)" = "phe",
                    "Tyrosine (tyr)" = "tyr",
                    "Acetate (ace)" = "ace",
                    "Acetoacetate (acace)" = "acace",
                    "Creatinine (crea)" = "crea",
                    "Albumin (alb)" = "alb",
                    "Glycoprotein acetyls (gp)" = "gp",
                    "Total lipids in chylomicrons and extremely large VLDL (xxl_vldl_l)" = "xxl_vldl_l",
                    "Cholesteryl esters in chylomicrons and extremely large VLDL (l_vldl_ce_percentage)"="l_vldl_ce_percentage",
                    "Free cholesterol in large LDL (l_hdl_fc)"="l_hdl_fc", 
                    "Ratio of apolipoprotein B to apolipoprotein A1 (apob_apoa1)"="apob_apoa1",
                    "Ratio of omega-6 fatty acids to omega-3 fatty acids (faw6_faw3)"="faw6_faw3",
                    "Glycine"="glycine"
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
  ),
  # Tab Panel for the histograms of the metabolite distributions compared to BBMRI-nl
  tabPanel(
    title = "Metabolites compared to BBMRI-nl",
    HTML('<hr style="border-color: #0088cc;">'),
    sidebarLayout(
      position = "right",
      sidebarPanel(
        width = 3,
        #Radio button if you want to scale them or not
        radioButtons(
          inputId = "scale_met_BBMRI",
          label="Do you want to see the scores scaled?",
          choices=c("Not scaled" = FALSE,
                    "scaled" = TRUE
          ),
          selected = FALSE,
        ),
        #Selection of which metabolites to print
        selectInput(
          inputId = "metabo_BBMRI",
          label="Metabolites distributions:",
          choices=c("Total lipids in medium VLDL (m_vldl_l)" = "m_vldl_l",
                    "Total lipids in small VLDL (s_vldl_l)" = "s_vldl_l",
                    "Total lipids in very small VLDL (xs_vldl_l)" = "xs_vldl_l",
                    "Total lipids in IDL (idl_l)" = "idl_l",
                    "Total cholesterol in IDL (idl_c)" = "idl_c",
                    "Total lipids in large LDL (l_ldl_l)" = "l_ldl_l",
                    "Total lipids in medium LDL (m_ldl_l)" = "m_ldl_l",
                    "Total lipids in medium HDL (m_hdl_l)" = "m_hdl_l",
                    "Total lipids in small HDL (s_hdl_l)" = "s_hdl_l",
                    "Total lipids in small LDL (s_ldl_l)" = "s_ldl_l",
                    "Mean diameter for VLDL particles (vldl_d)" = "vldl_d",
                    "Mean diameter for LDL particles (ldl_d)" = "ldl_d",
                    "Mean diameter for HDL particles (hdl_d)" = "hdl_d",
                    "Serum total cholesterol (serum_c)" = "serum_c",
                    "Total cholesterol in VLDL (vldl_c)" = "vldl_c",
                    "Total cholesterol in LDL (ldl_c)" = "ldl_c",
                    "Total cholesterol in HDL (hdl_c)" = "hdl_c",
                    "Total cholesterol in HDL2 (hdl2_c)" = "hdl2_c",
                    "Total cholesterol in HDL3 (hdl3_c)" = "hdl3_c",
                    "Total phosphoglycerides (totpg)" = "totpg",
                    "Phosphatidylcholine and other cholines (pc)" = "pc",
                    "Sphingomyelins (sm)" = "sm",
                    "Total cholines (totcho)" = "totcho",
                    "Apolipoprotein A-I (apoa1)" = "apoa1",
                    "Apolipoprotein B (apob)" = "apob",
                    "Total fatty acids (totfa)" = "totfa",
                    "Estimated degree of unsaturation (unsatdeg)" = "unsatdeg",
                    "docosahexaenoic acid (dha)" = "dha",
                    "linoleic acid (la)" = "la",
                    "Omega-3 fatty acids (faw3)" = "faw3",
                    "Omega-6 fatty acids (faw6)" = "faw6",
                    "Polyunsaturated fatty acids (pufa)" = "pufa",
                    "Monounsaturated fatty acids (mufa)" = "mufa",
                    "Saturated fatty acids (sfa)" = "sfa",
                    "Ratio of omega-3 fatty acids to total fatty acids (faw3_fa)" = "faw3_fa",
                    "Ratio of omega-6 fatty acids to total fatty acids (totfa)" = "totfa",
                    "Ratio of polyunsaturated fatty acids to total fatty acids (pufa_fa)" = "pufa_fa",
                    "Ratio of monounsaturated fatty acids to total fatty acids (mufa_fa)" = "mufa_fa",
                    "Ratio of saturated fatty acids to total fatty acids (sfa_fa)" = "sfa_fa",
                    "Glucose (glc)" = "glc",
                    "Lactate (lac)" = "lac",
                    "Citrate (cit)" = "cit",
                    "Alanine (ala)" = "ala",
                    "Glutamine (gln)" = "gln",
                    "Histidine (his)" = "his",
                    "Isoleucine (ile)" = "ile",
                    "Leucine (leu)" = "leu",
                    "Valine (val)" = "val",
                    "Phenylalanine (phe)" = "phe",
                    "Tyrosine (tyr)" = "tyr",
                    "Acetate (ace)" = "ace",
                    "Acetoacetate (acace)" = "acace",
                    "Creatinine (crea)" = "crea",
                    "Albumin (alb)" = "alb",
                    "Glycoprotein acetyls (gp)" = "gp",
                    "Total lipids in chylomicrons and extremely large VLDL (xxl_vldl_l)" = "xxl_vldl_l",
                    "Cholesteryl esters in chylomicrons and extremely large VLDL (l_vldl_ce_percentage)"="l_vldl_ce_percentage",
                    "Free cholesterol in large LDL (l_hdl_fc)"="l_hdl_fc", 
                    "Ratio of apolipoprotein B to apolipoprotein A1 (apob_apoa1)"="apob_apoa1",
                    "Ratio of omega-6 fatty acids to omega-3 fatty acids (faw6_faw3)"="faw6_faw3",
                    "Glycine"="glycine"
          ),
          selected = "m_vldl_l",
        ),
        style = "text-align: left;"
      ),
    #   mainPanel(
    #     width = 9,
    #     plotlyOutput("hist_BBMRI", height = "500px") %>% withSpinner()
    #   )
    # ),
    mainPanel(
      width = 9,
      fluidPage(
        fluidRow(
          #Histogram of the calibrations
          plotlyOutput("hist_BBMRI", height = "400px") %>% withSpinner()
        ),
        fluidRow(
          #Histogram of the calibrations
          plotlyOutput("models_coef_heat", height = "500px") %>% withSpinner()
        )
      )
    )
  ),
    
    helpText("Histograms to show the distributions of the metabolites uploaded compared to its distribution in BBMRI-nl.
             You can chose only one metabolite at a time. 
             You might want to look them scaled together in this comparison."),
    HTML('<hr style="border-color: #0088cc;">'),
  )
)
)
