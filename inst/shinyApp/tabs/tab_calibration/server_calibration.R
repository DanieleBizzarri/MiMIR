# Creates bot the calibration plots if the phenotype is available
calibration_plot<-reactive({
  bin<-as.character(input$surro_cal)
  sur<-paste0("s_",bin)
  ind<-which(MiMIR::phenotypes_names$bin_surro == sur)
  
  surro<-calculate_surrogate_scores(met=metabo_measures(), PARAM_surrogates = MiMIR::PARAM_surrogates,
                                    Nmax_miss=input$Nmax_miss_surrogates,
                                    Nmax_zero=input$Nmax_zero_surrogates,
                                    bin_names = MiMIR::phenotypes_names$bin_names,
                                    roc=F, quiet=T, post=F)
  surro<-surro$surrogates
  if(is.null(calibrations()[[ind]])){
    return(NULL)
  }else{
    return(
      suppressWarnings(plattCalib_evaluation(r=bin_phenotypes()[,bin], p=surro[,sur],p.orig=surrogates()[,sur],
                            name=sur, nbins = input$Nbins)))
  }
})

# Ouptuts reliability plot if the phenotype is available
output$reliability_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      if(is.null(calibration_plot())){
        suppressWarnings(plotly_NA_message(main="This phenotype is not available"))
      }else{
        suppressWarnings(calibration_plot()$cal.Plot)
      }
      }, error = function(err) {
        return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
      })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Ouptuts histogram of the calibrations if the phenotype is available
output$hist_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
    if(is.null(calibration_plot())){
      plotly_NA_message(main="This phenotype is not available")
    }else{
      suppressWarnings(calibration_plot()$prob.hist)
    }
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Ouptuts the heatmap for the correlations of the calibrated surrogates
output$heat_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(calibrations())
      cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
      res<-cor_assoc(cal,cal, colnames(cal), colnames(cal))
      heat<-plot_corply(res, main="Calibrated surrogates Correlations", reorder.x=TRUE, abs=F, 
                        resort_on_p= TRUE,reorder_dend=F)
      
    }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

# Ouptuts the missingness heatmap for the calibrated surrogates
output$heat_na_calib <- renderPlot({
  if(required()){
    tryCatch({
      req(calibrations())
      cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
      plot_na_heatmap(t(cal))
    }, error = function(err) {
      return(NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})