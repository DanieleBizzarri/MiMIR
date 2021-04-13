output$hist_predictors <- renderPlotly({
  tryCatch({
    req(required())
    pred<-as.character(input$predictors)
    hist_plots(predictors(),x_name=pred, scaled=input$scaling)
  }, error = function(err) {
    return(met_NA())
  })
})

output$heat_predictors <- renderPlotly({
  tryCatch({
    res<-cor.assoc(predictors(),predictors(), colnames(predictors())[-1],colnames(predictors())[-1])
    heat<-plot.corply(res, main="Correlation of the Metabolic scores", reorder.x=TRUE, abs=F, 
                      resort_on_p= TRUE,reorder_dend=F)
  }, error = function(err) {
    return(met_NA())
  })
})

output$heat_na_pred <- renderPlot({
  if(required()){
      plot.na.heatmap(t(predictors()))
  }else{
    return(met_NA())
  }
})

