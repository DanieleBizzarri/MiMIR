#Creates the plotly heatmap with the phenotypes' correlations
output$heat_pheno <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      res<-cor_assoc(data.matrix(bin_phenotypes()[,bin_pheno_available()]),data.matrix(bin_phenotypes()[,bin_pheno_available()]), bin_pheno_available(), bin_pheno_available())
      heat<-plot_corply(res, main="Binary Phenotypes Correlations", reorder.x=TRUE, abs=F, 
                        resort_on_p= TRUE,reorder_dend=F)
      
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

#Creates the figure with the phenotypes' missingness
output$heat_na_pheno <- renderPlot({
  if(required()){
    tryCatch({
      req(phenotypes())
      plot_na_heatmap(t(bin_phenotypes()))
    }, error = function(err) {
      return(NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

output$summ_pheno <-  DT::renderDataTable({
  tryCatch({
    req(phenotypes())
    summary(phenotypes())
  }, error = function(err) {
    return(NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
  })
})

output$bin_phen_barplot <- renderPlotly({
  tryCatch({
    req(phenotypes())
    pheno_barplots(bin_phenotypes())
  }, error = function(err) {
    return(NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
  })
})


