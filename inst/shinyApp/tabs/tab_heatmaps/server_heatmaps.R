# output$heat_met <- renderPlotly({
#   tryCatch({
#     req(required())
#   if (input$MET56_14_cor == "MET56") {
#     met = MET56
#   } else {
#     met = MET14
#   }
#   res<-cor.assoc(metabo_measures(),metabo_measures(), met,met)
#   heat<-plot.corply(res, main="Metabolites' Correlations", reorder.x=TRUE, abs=F, 
#                     resort_on_p= TRUE,reorder_dend=F)
#   return(heat)
#   }, error = function(err) {
#     return(met_NA())
#   })
# })
# 
# output$heat_na_metabo <- renderPlot({
#   tryCatch({
#     req(required())
#     if (input$MET56_14_nas == "MET56") {
#       met = MET56
#     } else {
#       met = MET14
#     }
#     plot.na.heatmap(t(metabo_measures()[,met]))
#   }, error = function(err) {
#     return(met_NA_image())
#   })
# })


# output$heat_na_pheno <- renderPlot({
#   if(required()){
#     tryCatch({
#       req(phenotypes())
#     plot.na.heatmap(t(bin_phenotypes()))
#     }, error = function(err) {
#       return(pheno_NA_image())
#     })
#   }else{
#     return(met_NA())
#   }
# })

output$heat_predictors <- renderPlotly({
  tryCatch({
  res<-cor.assoc(predictors(),predictors(), colnames(predictors())[-1],colnames(predictors())[-1])
  heat<-plot.corply(res, main="Correlation of the Metabolic scores", reorder.x=TRUE, abs=F, 
                    resort_on_p= TRUE,reorder_dend=F)
  }, error = function(err) {
    return(met_NA())
  })
})
