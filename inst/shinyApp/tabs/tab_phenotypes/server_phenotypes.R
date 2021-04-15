#Creates the plotly heatmap with the phenotypes' correlations
output$heat_pheno <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      res<-cor.assoc(data.matrix(bin_phenotypes()[,bin_pheno_available()]),data.matrix(bin_phenotypes()[,bin_pheno_available()]), bin_pheno_available(), bin_pheno_available())
      heat<-plot.corply(res, main="Binary Phenotypes' Correlations", reorder.x=TRUE, abs=F, 
                        resort_on_p= TRUE,reorder_dend=F)
      
    }, error = function(err) {
      return(pheno_NA_image())
    })
  }else{
    return(met_NA())
  }
})

#Creates the figure with the phenotypes' missingness
output$heat_na_pheno <- renderPlot({
  if(required()){
    tryCatch({
      req(phenotypes())
      plot.na.heatmap(t(bin_phenotypes()))
    }, error = function(err) {
      return(pheno_NA_image())
    })
  }else{
    return(met_NA())
  }
})

