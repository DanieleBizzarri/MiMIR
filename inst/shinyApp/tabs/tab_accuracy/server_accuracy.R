# Creates the ROC curves for the selected surrogates
output$ROCs <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      surrog<-as.character(input$surroc)
      roc_surro(surrogates(),bin_phenotypes(),x_name=surrog)
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
    }else{
      return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Show the t-tests for the surrogates in the dataset
output$ttest <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      #ttest_surrogates(surrogates = predictors()[,bin_surro], bin_phenotypes = bin_phenotypes())
      ttest_surrogates(surrogates = surrogates(), bin_phenotypes = bin_phenotypes())
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Show the t-tests for the surrogates in the dataset
output$LOBOV_surro <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      LOBOV_accuracies(surrogates= surrogates(), bin_phenotypes= bin_phenotypes(), bin_pheno_available = bin_pheno_available(), acc_LOBOV= MiMIR::acc_LOBOV)
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

