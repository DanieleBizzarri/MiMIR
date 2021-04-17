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
      ttest_surrogates(surrogates = surrogates(), bin_phenotypes = bin_phenotypes())
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Creates a scatterplot comparing MetaboAge and age
output$scatter_metaboage <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      if(length(which(colnames(phenotypes())=="age"))==0){
        plotly_NA_message(main="Age is not available")
      }else{
        x<-data.frame(phenotypes()[,"age"])
        rownames(x)<-rownames(phenotypes())
        scatterplot_predictions(x, MetaboAge(), 
                                title="Chronological Age vs MetaboAge",
                                xname="Chronological age", 
                                yname="MetaboAge")
        }
      }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

#Creates an histogram of the mortality score divided for the different ages
output$hist_mort <- renderPlotly({
  if(required()){
      if(dim(phenotypes())[1]!=0){
        if(length(which(colnames(phenotypes())=="age"))==0){
          plotly_NA_message(main="age is not available")
        }else{
          hist_plots_mortality(mort_score(),phenotypes())
        }
      }else{
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    }
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})
