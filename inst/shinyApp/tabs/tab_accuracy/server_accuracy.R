output$ROCs <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      surrog<-as.character(input$surroc)
      roc_surro(surrogates(),bin_phenotypes(),x_name=surrog)
    }, error = function(err) {
      return(phenos_NA())
    })
    }else{
      return(met_NA())
  }
})

output$ttest <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      ttest_surrogates(surrogates = surrogates(), bin_phenotypes = bin_phenotypes())
    }, error = function(err) {
      return(phenos_NA())
    })
  }else{
    return(met_NA())
  }
})


output$scatter_metaboage <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      if(length(which(colnames(phenotypes())=="age"))==0){
        pheno_NA()
      }else{
        x<-data.frame(phenotypes()[,"age"])
        rownames(x)<-rownames(phenotypes())
        scatterplot_predictions(x, MetaboAge(), 
                                title="Chronological Age vs MetaboAge",
                                xname="Chronological age", 
                                yname="MetaboAge")
        }
      }, error = function(err) {
      return(phenos_NA())
    })
  }else{
    return(met_NA())
  }
})


output$hist_mort <- renderPlotly({
  if(required()){
      if(dim(phenotypes())[1]!=0){
        if(length(which(colnames(phenotypes())=="age"))==0){
          pheno_NA()
        }else{
          hist_plots_mortality(mort_score(),phenotypes())
        }
      }else{
      return(phenos_NA())
    }
  }else{
    return(met_NA())
  }
})
