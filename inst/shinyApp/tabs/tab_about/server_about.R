## Description of our works and how to cite them
output[["description"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "This R-shiny webtool was built to be able to compute the multi-variate metabolomisc scores from raw Nightingale Health 1H-NMR metabolomics data.
        Please refer to the original manuscripts when using these metabolic biomarkers in your works:
        ",
        style = "text-align: justify; display: block; width: 90%"
      ),
      span(
        "1) mortality score: Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Commu-nications, 10, 1–8.",
        tags$br(),
        "2) MetaboAge: van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547.",
        tags$br(),
       "3) surrogate clinical variables: Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764.",
       tags$br(),
       "4) COVID score: Nightingale Health UK Biobank Initiative et al. (2021) Metabolic biomarker profiling for identification of susceptibility to severe pneumonia and COVID-19 in the general population. eLife, 10, e63033",
       tags$br(),
       "5) T2D score: Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298–2309.",
       tags$br(),
       "6) CVD score: Würtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774–785.",
       style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
    )
  )
})


## Get sessionInfo from current session
output[["currentSession"]] <- renderPrint({
  current <- sessionInfo()
  current
})