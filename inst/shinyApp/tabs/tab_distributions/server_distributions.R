# Outputs the histograms of the predictors chosen (scaled or not)
output$hist_predictors <- renderPlotly({
  tryCatch({
    req(required())
    pred<-as.character(input$predictors)
    hist_plots(predictors(),x_name=pred, scaled=input$scaling)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

# Outputs the correlation heatmap of the predictors
output$heat_predictors <- renderPlotly({
  tryCatch({
    res<-cor_assoc(predictors(),predictors(), colnames(predictors())[-1],colnames(predictors())[-1])
    heat<-plot_corply(res, main="Correlation of the Metabolic scores", reorder.x=TRUE, abs=F, 
                      resort_on_p= TRUE,reorder_dend=F)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

# Outputs the heatmap of the missingness in the predictors
output$heat_na_pred <- renderPlot({
  if(required()){
      plot_na_heatmap(t(predictors()))
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

