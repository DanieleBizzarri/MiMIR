## QC introduction 
output[["qc_intro"]] <- renderUI({
  tagList(
    div(
      span(
        "The mortality score, MetaboAge and the surrogate clinical values were built in different moments
        and with different strategies, therefore they have a different Quality Control Protocol.
        In this page we go through the QC for each of the methods.
        ",
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      span(
        paste("The mortality score has a fixed QC, which  consists in the log transform of the z-scaled metabolites done separately for each study.
        It doesn't involve  any imputation, therefore the missingness in your dataset relative to the 14 metabolites present in the score will cause
        the missingness of the mostality score value. The number of available mortality score values in the current dataset is:
        ", length(which(!is.na(comp.mort_score(metabo_measures(), quiet=TRUE))))),
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      span(
        "For the MetaboAge and the surrogates scores instead we z-scale the metabolites and impute with zero the missing values.
        We allow the samples that don't have metabolites measures going over 5 times the standard deviation we measured in BBMRI-NL,
        and we allow by default: 1 missing value per sample (Nmax_miss) and 1 zeros per sample (Nmax_zero).
        We reccomend to use this values, but feel free to change the values of Nmax_miss and Nmax_zero separately for MetaboAge and the surrogates.
        
        ",
        style = "text-align: justify; display: block; width: 90%"
      ),
      span(
        "When you decided the best for you, press the button \"Run Analysis!\" at the bottom of the page!",
        style = "text-align: justify; display: block; width: 90%"
      )
    )
  )
})

# Output MetaboAge settings
output$MetaboAge_settings <- renderUI({
  str1 <- "You have selected the current values:"
  str2 <- paste0("N max zeros= ", input$Nmax_zero_metaboAge,";")
  str3 <- paste0("N max missing= ", input$Nmax_miss_metaboAge,".")
  
  HTML(paste(str1, str2, str3, sep = '<br/>'))
})

# Output Surrogates settings
output$Surrogates_settings <- renderUI({
  str1 <- "You have selected the current values:"
  str2 <- paste0("N max zeros= ", input$Nmax_zero_surrogates,";")
  str3 <- paste0("N max missing= ", input$Nmax_miss_surrogates,".")
  
  HTML(paste(str1, str2, str3, sep = '<br/>'))
})

# Output MetaboAge resulting dataset info
output$QC_metaboAge_text <- renderPrint(result <- QCprep(as.matrix(metabo_measures()[,MET63]),
                                                         PARAM,quiet=FALSE,
                                                         Nmax_miss=input$Nmax_miss_metaboAge,
                                                         Nmax_zero=input$Nmax_zero_metaboAge))
# Output Surrogates resulting dataset info
output$QC_surrogates_text <- renderPrint(result <- QCprep_surrogates(as.matrix(metabo_measures()[,MET63]),
                                                                      PARAM_surrogates,quiet=FALSE,
                                                                      Nmax_miss=input$Nmax_miss_surrogates,
                                                                      Nmax_zero=input$Nmax_zero_surrogates))

# Calculate the mortality score
mort_score <- reactive({
  return(comp.mort_score(metabo_measures(),quiet=TRUE))
})

# Calculate the MetaboAge
MetaboAge <- reactive({
  metabo_metaboage<-QCprep(as.matrix(metabo_measures()[,MET63]),
                           PARAM,quiet=TRUE,
                           Nmax_miss=input$Nmax_miss_metaboAge,
                           Nmax_zero=input$Nmax_zero_metaboAge)
  
  return(apply.fit(metabo_metaboage,FIT=PARAM$FIT_COEF))
})

# Calculate the surrogates
surrogates <- reactive({
  surro<-calculate_surrogate_scores(met=metabo_measures(), PARAM_surrogates = PARAM_surrogates,
                                    Nmax_miss=input$Nmax_miss_surrogates,
                                    Nmax_zero=input$Nmax_zero_surrogates,
                                    bin_names = bin_names,
                                    roc=F, quiet=T)
  return(surro$surrogates)
})

# Compose the predictors table with mortality score, metaboage and surrogates
predictors<-reactive({
  mortality<-data.frame(ID=row.names(mort_score()),mort_score())
  metaboage<-data.frame(ID=rownames(MetaboAge()),MetaboAge())
  surro<-data.frame(ID=rownames(surrogates()),surrogates())
  
  mort_age<-merge(mortality,metaboage, by.x= 'ID',by.y="ID", all=TRUE)
  predictors<-merge(mort_age,surro, by.x= 'ID',by.y="ID", all=TRUE)
  return(predictors)
})

# Calibrated surrogates
calibrations<-reactive({
  calib<-calibration_surro(bin_phenotypes=bin_phenotypes(), surrogates=surrogates(), 
                           bin_names=bin_names, bin_pheno_available=bin_pheno_available(), pl=FALSE)
})

