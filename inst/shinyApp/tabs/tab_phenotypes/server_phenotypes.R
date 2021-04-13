output$heat_pheno <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      #res<-cor.assoc(data.matrix(phenotypes()[,pheno_available()]),data.matrix(phenotypes()[,pheno_available()]), pheno_available(), pheno_available())
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


# output$hist_pheno <- renderPlotly({
#   if(required()){
#     tryCatch({
#       req(phenotypes())
#       pheno<-as.character(input$pheno)
#       hist_plots(bin_phenotypes(),x_name=pheno, scaled=F)
#     }, error = function(err) {
#       return(pheno_NA_image())
#     })
#   }else{
#     return(met_NA())
#   }
# })