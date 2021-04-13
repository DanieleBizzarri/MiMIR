## QC introduction test
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


# source("script/predictors_functions.R", local = TRUE)
# load("data/PARAM__2018-06-18_02-16-17.457.RData")
# load("data/PARAM_surrogates_2021_02_23.RData")

output$MetaboAge_settings <- renderUI({
  str1 <- "You have selected the current values:"
  str2 <- paste0("N max zeros= ", input$Nmax_zero_metaboAge,";")
  str3 <- paste0("N max missing= ", input$Nmax_miss_metaboAge,".")
  
  HTML(paste(str1, str2, str3, sep = '<br/>'))
})

output$Surrogates_settings <- renderUI({
  str1 <- "You have selected the current values:"
  str2 <- paste0("N max zeros= ", input$Nmax_zero_surrogates,";")
  str3 <- paste0("N max missing= ", input$Nmax_miss_surrogates,".")
  
  HTML(paste(str1, str2, str3, sep = '<br/>'))
})


output$QC_metaboAge_text <- renderPrint(result <- QCprep(as.matrix(metabo_measures()[,MET63]),
                                                         PARAM,quiet=FALSE,
                                                         Nmax_miss=input$Nmax_miss_metaboAge,
                                                         Nmax_zero=input$Nmax_zero_metaboAge))

output$QC_surrogates_text <- renderPrint(result <- QCprep_metabotypes(as.matrix(metabo_measures()[,MET63]),
                                                                      PARAM_surrogates,quiet=FALSE,
                                                                      Nmax_miss=input$Nmax_miss_surrogates,
                                                                      Nmax_zero=input$Nmax_zero_surrogates))


mort_score <- reactive({
  return(comp.mort_score(metabo_measures(),quiet=TRUE))
})

MetaboAge <- reactive({
  metabo_metaboage<-QCprep(as.matrix(metabo_measures()[,MET63]),
                           PARAM,quiet=TRUE,
                           Nmax_miss=input$Nmax_miss_metaboAge,
                           Nmax_zero=input$Nmax_zero_metaboAge)
  
  return(apply.fit(metabo_metaboage,FIT=PARAM$FIT_COEF))
})

surrogates <- reactive({
  surro<-calculate_surrogate_scores(met=metabo_measures(), PARAM_surrogates = PARAM_surrogates,
                                    Nmax_miss=input$Nmax_miss_surrogates,
                                    Nmax_zero=input$Nmax_zero_surrogates,
                                    bin_names = bin_names,
                                    roc=F, quiet=T)
  return(surro$surrogates)
})


predictors<-reactive({
  mortality<-data.frame(ID=row.names(mort_score()),mort_score())
  metaboage<-data.frame(ID=rownames(MetaboAge()),MetaboAge())
  surro<-data.frame(ID=rownames(surrogates()),surrogates())
  
  a<-merge(mortality,metaboage, by.x= 'ID',by.y="ID", all=TRUE)
  predictors<-merge(a,surro, by.x= 'ID',by.y="ID", all=TRUE)
  return(predictors)
})

calibrations<-reactive({
  # calib<-lapply(1:length(bin_surro), function(i){
  #   orig<-as.numeric(bin_phenotypes()[,bin_names[i]])-1
  #   pred<-as.numeric(surrogates()[,bin_surro[i]])
  #   ind<-which(!is.na(orig))
  #   if(length(ind)!=0){
  #     calibration<-plattCalibration(orig[ind],pred[ind], 10)
  #     return(calibration)
  #   }
  # })
  calib<-calibration_surro(bin_phenotypes(), surrogates(), bin_names, bin_surro, bin_pheno_available())
})

# ## Start an analysis
# observeEvent(input$run_button, {
#   
#   mort_score <- reactive({
#     comp.mort_score(metabo_measures())
#   })
#   
#  MetaboAge <- reactive({
#    metabo_metaboage<-QCprep(as.matrix(metabo_measures()[,MET63]),
#           PARAM,quiet=TRUE,
#           Nmax_miss=input$Nmax_miss_metaboAge,
#           Nmax_zero=input$Nmax_zero_metaboAge)
#    
#    return(apply.fit(metabo_metaboage,FIT=PARAM$FIT_COEF))
#   })
#   
#   
# })

