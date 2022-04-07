## QC introduction 
output[["qc_intro"]] <- renderUI({
  if (is.null(syn_met()) & is.null(upload_met())) {
    tagList(div(span("Please first upload the Nightingale metabolomics features!",
                     style = "text-align: justify; display: block; width: 90%")))
    }else{
    tagList(
      div(
        span(
          "The metabolomics-based scores were built using different pre-processing strategies. 
          The Quality Control steps are, for the most part already fixed by the original authors of the scores. 
          This process generally consists in scaling the metabolic features, to selecting samples that were correctly measured and eventually impute the missing values.
          This tab gives the option to make the user able to have freedom, if possible, over some Steps in some metabolomics-bases scores.
          ",
          style = "text-align: justify; display: block; width: 90%"
        ),
        br(),
        span(
          paste("Mortality, COVID, T2D and CVD scores have a fixed QC, which  consists in log transforming the z-scaled metabolites.
          It doesn't involve  any imputation, therefore the missingness in your dataset will cause missingness in the scores.
          "),
          style = "text-align: justify; display: block; width: 90%"
        ),
        br(),
        span(
          "MetaboAge and the surrogates scores: We select samples with limited number of zeros, missing values and exclude samples with outliers (based on the means and standard deviations in BBMRI-NL),
          MiMIR calculates the scores automatically with maximum 1 missing value (Nmax_miss), 1 zeros (Nmax_zero) per sample.
          The user can to change the values of Nmax_miss and Nmax_zero accordingly to his/her needs.
          Finally, the metabolites are z-scaled and imputed with zero the missing values.
          ",
          style = "text-align: justify; display: block; width: 90%"
        ),
        br(),
        span(
          "Please look at the MANUAL for a complete description of the Pre-processing steps in all the scores.
          ",
          style = "text-align: justify; display: block; width: 90%"
        )
      )
    )}
  })
  
# Output MetaboAge settings
output$MetaboAge_settings <- renderUI({
  if(required()){
      str1 <- "You have selected the current values:"
      str2 <- paste0("N max zeros= ", input$Nmax_zero_metaboAge,";")
      str3 <- paste0("N max missing= ", input$Nmax_miss_metaboAge,".")
      HTML(paste(str1, str2, str3, sep = '<br/>'))
}else{
  div(HTML("Please Upload the Nightingale metabolomics features first!"))
}
})

# Output Surrogates settings
output$Surrogates_settings <- renderUI({
  if(required()){
  str1 <- "You have selected the current values:"
  str2 <- paste0("N max zeros= ", input$Nmax_zero_surrogates,";")
  str3 <- paste0("N max missing= ", input$Nmax_miss_surrogates,".")
  
  HTML(paste(str1, str2, str3, sep = '<br/>'))
  }else{
    tagList(div(span("Please Upload the Nightingale metabolomics features first!",
                     style = "text-align: justify; display: block; width: 90%")))
  }
})

# Output MetaboAge resulting dataset info
output$QC_metaboAge_text <- renderPrint(
  if(required()){
  result <- QCprep(as.matrix(metabo_measures()[,MiMIR::metabolites_subsets$MET63]),
                                                         MiMIR::PARAM_metaboAge,quiet=FALSE,
                                                         Nmax_miss=input$Nmax_miss_metaboAge,
                                                         Nmax_zero=input$Nmax_zero_metaboAge)
  }else{
    return(cat("Please Upload the Nightingale metabolomics features first."))
  })

# Output Surrogates resulting dataset info
output$QC_surrogates_text <- renderPrint(
  if(required()){
  result <- QCprep_surrogates(as.matrix(metabo_measures()[,MiMIR::metabolites_subsets$MET63]),
                                                                     MiMIR::PARAM_surrogates,quiet=FALSE,
                                                                     Nmax_miss=input$Nmax_miss_surrogates,
                                                                     Nmax_zero=input$Nmax_zero_surrogates)
  }else{
    return("Please Upload the Nightingale metabolomics features first.")
  })

# Calculate the mortality score
mort_score <- reactive({
  mortsc<-comp.mort_score(metabo_measures(),quiet=TRUE)
  return(mortsc)
})

# Calculate the MetaboAge
MetaboAge <- reactive({
  metabo_metaboage<-QCprep(as.matrix(metabo_measures()[,MiMIR::metabolites_subsets$MET63]),
                           MiMIR::PARAM_metaboAge,quiet=TRUE,
                           Nmax_miss=input$Nmax_miss_metaboAge,
                           Nmax_zero=input$Nmax_zero_metaboAge)
  metaboage<-apply.fit(metabo_metaboage,FIT=MiMIR::PARAM_metaboAge$FIT_COEF)
  
  #rownames(metaboage)<-rownames(metabo_measures())
  if(!is.null(phenotypes()[,"age"])){
    metaboage<-data.frame(metaboage=metaboage, deltaMetaboAge=metaboage-(phenotypes()[,"age"]))
    colnames(metaboage)<-c("MetaboAge","deltaMetaboAge")
    return(metaboage)
  }else{
    return(metaboage)
  }
})

# Calculate the surrogates
surrogates <- reactive({
  surro<-calculate_surrogate_scores(met=metabo_measures(), 
                                    PARAM_surrogates = MiMIR::PARAM_surrogates,
                                    Nmax_miss=input$Nmax_miss_surrogates,
                                    Nmax_zero=input$Nmax_zero_surrogates,
                                    bin_names = MiMIR::phenotypes_names$bin_names,
                                    roc=F, quiet=T)
  surrogates<-surro$surrogates
  return(surrogates)
})

# Calculate the T2D score
T2D_score_AholaOlli <- reactive({
  T2Dscore<-comp.T2D_Ahola_Olli(met=metabo_measures(), phen=phenotypes(),quiet=TRUE)
})
# Calculate CVD score
CVD_score <- reactive({
  CVDscore<-comp.CVD_score(met=metabo_measures(), phen=phenotypes(),quiet=TRUE)
})

# Calculate CVD score
COVID_score <- reactive({
  COVIDscore<-comp_covid_score(dat=metabo_measures(), quiet=TRUE)
})

# Compose the predictors table with mortality score, metaboage and surrogates
predictors<-reactive({
  mortality<-data.frame(ID=rownames(mort_score()), mort_score())
  metaboage<-data.frame(ID=rownames(MetaboAge()), MetaboAge())
  surro<-data.frame(ID=rownames(surrogates()), surrogates())
  
  predictors<-merge(mortality,metaboage, by.x= 'ID',by.y="ID", all=TRUE)
  if(length(t(T2D_score_AholaOlli()))>2){
    T2D_score<-data.frame(ID=rownames(T2D_score_AholaOlli()),T2D_score_AholaOlli())
    predictors<-merge(predictors,T2D_score, by.x= 'ID',by.y="ID", all=TRUE)
  }
  if(length(t(CVD_score()))>2){
    CVDscore<-data.frame(ID=rownames(CVD_score()),CVD_score())
    predictors<-merge(predictors,CVDscore, by.x= 'ID',by.y="ID", all=TRUE)
  }
  if(length(t(COVID_score()))>2){
    COVID_score<-data.frame(ID=rownames(COVID_score()),COVID_score())
    predictors<-merge(predictors,COVID_score, by.x= 'ID',by.y="ID", all=TRUE)
  }
  
  predictors<-merge(predictors,surro, by.x= 'ID',by.y="ID", all=TRUE)
  
  ind<-order(match(predictors$ID, rownames(metabo_measures())))
  predictors<-predictors[ind,]
  rownames(predictors)<-predictors$ID
  
  return(predictors)
})

# Calibrated surrogates
calibrations<-reactive({
  surro<-calculate_surrogate_scores(met=metabo_measures(), 
                                    PARAM_surrogates = MiMIR::PARAM_surrogates,
                                    Nmax_miss=input$Nmax_miss_surrogates,
                                    Nmax_zero=input$Nmax_zero_surrogates,
                                    bin_names = MiMIR::phenotypes_names$bin_names,
                                    roc=F, quiet=T, post=F)
  surro<-surro$surrogates
  calib<-calibration_surro(bin_phenotypes=bin_phenotypes(), surrogates=surro, 
                           bin_names= MiMIR::phenotypes_names$bin_names, bin_pheno_available=bin_pheno_available(), pl=FALSE)
})

