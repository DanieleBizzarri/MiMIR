#Creates the plotly heatmap with the metabolites' correlations
output$heat_met <- renderPlotly({
  tryCatch({
    req(required())
    if (input$MET56_14_cor == "MET56") {
      met = MET56
    } else if (input$MET56_14_cor == "MET14"){
      met = MET14
    } else if (input$MET56_14_cor == "MET_COVID"){
      met = MET_COVID
    }else if (input$MET56_14_cor == "MET_T2D"){
      met = MET_T2D
    }else if (input$MET56_14_cor == "MET_CVD"){
      met = MET_CVD
    }
    
    if(any(colSums(is.na(metabo_measures()[,met])) %in% dim(metabo_measures())[1])){
      met<-met[-which(colSums(is.na(metabo_measures()[,met])) %in% dim(metabo_measures())[1])]
    }
    
    res<-cor_assoc(metabo_measures(),metabo_measures(), met, met)
    heat<-plot_corply(res, main="Metabolites' Correlations", reorder.x=TRUE, abs=F, 
                      resort_on_p= TRUE,reorder_dend=F)
    return(heat)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

#Creates the figure with the metabolites' missingness
output$heat_na_metabo <- renderPlot({
  tryCatch({
    req(required())
    if (input$MET56_14_nas == "MET56") {
      met = MET56
    } else if (input$MET56_14_nas == "MET14"){
      met = MET14
    } else if (input$MET56_14_nas == "MET_COVID"){
      met = MET_COVID
    }else if (input$MET56_14_nas == "MET_T2D"){
      met = MET_T2D
    }else if (input$MET56_14_nas == "MET_CVD"){
      met = MET_CVD
    }
    
    plot_na_heatmap(t(metabo_measures()[,met]))
  }, error = function(err) {
    return(NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

#Creates the plotly histogram of the selected metabolites
output$hist_metabolites <- renderPlotly({
  tryCatch({
    req(required())
    met<-as.character(input$metabolites)
    hist_plots(metabo_measures(),x_name=met, scaled=input$scale_met, 
               main="Metabolites Distributions",datatype="metabolites")
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

#Creates the plotly histogram of the selected metabolites
output$hist_BBMRI <- renderPlotly({
  tryCatch({
    req(required())
    met<-as.character(input$metabo_BBMRI)
    BBMRI_hist_plot(dat=metabo_measures(), x_name=met, scaled=input$scale_met_BBMRI)
  }, error = function(err) {
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  })
})

#Creates the plotly histogram of the selected metabolites
output$models_coef_heat <- renderPlotly({
  model_coeff_heat(mort_betas=mort_betas,metaboAge_betas=PARAM$FIT_COEF, surrogates_betas=PARAM_surrogates$models_betas,
                  Ahola_Olli_betas= Ahola_Olli_betas, CVD_score_betas= CVD_score_betas, COVID_score_betas= covid_betas)
})

