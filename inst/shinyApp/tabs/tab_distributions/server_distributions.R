# Outputs the histograms of the surrogates chosen (scaled or not)
output$hist_surrogates <- renderPlotly({
  tryCatch({
    req(required())
    pred<-as.character(input$surrogates)
    hist_plots(surrogates(), x_name=pred, scaled=input$scaling_surro)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

# Outputs the correlation heatmap of the surrogates
output$heat_surrogates <- renderPlotly({
  tryCatch({
    res<-cor_assoc(surrogates(),surrogates(), colnames(surrogates())[-1],colnames(surrogates())[-1])
    heat<-plot_corply(res, main="Correlation of the Metabolic surrogates", reorder.x=TRUE, abs=F, 
                      resort_on_p= TRUE,reorder_dend=F)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

# Outputs the heatmap of the missingness in the surrogates
output$heat_na_surro <- renderPlot({
  if(required()){
      plot_na_heatmap(t(surrogates()))
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

