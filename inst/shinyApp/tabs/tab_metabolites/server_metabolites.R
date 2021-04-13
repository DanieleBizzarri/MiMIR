output$heat_met <- renderPlotly({
  tryCatch({
    req(required())
    if (input$MET56_14_cor == "MET56") {
      met = MET56
    } else {
      met = MET14
    }
    res<-cor.assoc(metabo_measures(),metabo_measures(), met,met)
    heat<-plot.corply(res, main="Metabolites' Correlations", reorder.x=TRUE, abs=F, 
                      resort_on_p= TRUE,reorder_dend=F)
    return(heat)
  }, error = function(err) {
    return(met_NA())
  })
})

output$heat_na_metabo <- renderPlot({
  tryCatch({
    req(required())
    if (input$MET56_14_nas == "MET56") {
      met = MET56
    } else {
      met = MET14
    }
    plot.na.heatmap(t(metabo_measures()[,met]))
  }, error = function(err) {
    return(met_NA_image())
  })
})

output$hist_metabolites <- renderPlotly({
  tryCatch({
    req(required())
    met<-as.character(input$metabolites)
    hist_plots(metabo_measures(),x_name=met, scaled=input$scale_met, main="Metabolites' Distributions")
  }, error = function(err) {
    return(met_NA())
  })
})