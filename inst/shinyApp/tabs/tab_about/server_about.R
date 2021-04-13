## Welcome text
output[["description"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "The R-shiny webtool to calculate the metabolic predictors from raw Nightingale Health 1H-NMR metabolomics data,
        developed by the groups of Molecular Epidemiology and LCBC (LUMC).
        Please refer to our manuscripts when using these metabolic biomarkers in your works:
        ",
        style = "text-align: justify; display: block; width: 90%"
      ),
      span(
        "1) mortality score: J. Deelen et al., ‘A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals’, Nat. Commun., vol. 10, no. 1, pp. 1–8, Aug. 2019, doi: 10.1038/s41467-019-11311-9",
        tags$br(),
        "2) MetaboAge: van den Akker Erik B. et al., ‘Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease’, Circ. Genomic Precis. Med., vol. 13, no. 5, pp. 541–547, Oct. 2020, doi: 10.1161/CIRCGEN.119.002610.",
        tags$br(),
       "3) surrogate clinical variables: unpublished",
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
    )
  )
})

output[["references"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
      "Please refer to our manuscripts when using these metabolic biomarkers in your works:",
      style = "text-align: justify; display: block; width: 90%"
      ),
      span(
        "[1] van den Akker Erik B. et al., ‘Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease’, Circ. Genomic Precis. Med., vol. 13, no. 5, pp. 541–547, Oct. 2020, doi: 10.1161/CIRCGEN.119.002610.",
        tags$br(),
        "[2] J. Deelen et al., ‘A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals’, Nat. Commun., vol. 10, no. 1, pp. 1–8, Aug. 2019, doi: 10.1038/s41467-019-11311-9",
        tags$br(),
        "[3] unpublished",
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