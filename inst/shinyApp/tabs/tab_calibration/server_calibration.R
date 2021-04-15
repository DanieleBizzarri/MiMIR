# Creates bot the calibration plots if the phenotype is available
calibration_plot<-reactive({
  sur<-as.character(input$surrogates)
  ind<-which(bin_surro == sur)
  if(is.null(calibrations()[[ind]])){
    return(NULL)
  }else{
    return(plattCalib_plot(calibrations()[[ind]],name=sur, nbins = input$Nbins))
  }
})

# Ouptuts reliability plot if the phenotype is available
output$reliability_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      if(is.null(calibration_plot())){
        pheno_NA()
      }else{
        calibration_plot()$cal.Plot
      }
      }, error = function(err) {
        return(phenos_NA())
      })
  }else{
    return(met_NA())
  }
})

# Ouptuts histogram of the calibrations if the phenotype is available
output$hist_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
    if(is.null(calibration_plot())){
      pheno_NA()
    }else{
      calibration_plot()$prob.hist
    }
    }, error = function(err) {
      return(phenos_NA())
    })
  }else{
    return(met_NA())
  }
})

# Ouptuts the heatmap for the correlations of the calibrated surrogates
output$heat_calib <- renderPlotly({
  if(required()){
    tryCatch({
      req(calibrations())
      cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
      res<-cor.assoc(cal,cal, colnames(cal), colnames(cal))
      heat<-plot.corply(res, main="Calibrated surrogates Correlations", reorder.x=TRUE, abs=F, 
                        resort_on_p= TRUE,reorder_dend=F)
      
    }, error = function(err) {
      return(pheno_NA_image())
    })
  }else{
    return(met_NA())
  }
})

# Ouptuts the missingness heatmap for the calibrated surrogates
output$heat_na_calib <- renderPlot({
  if(required()){
    tryCatch({
      req(calibrations())
      cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
      plot.na.heatmap(t(cal))
    }, error = function(err) {
      return(pheno_NA_image())
    })
  }else{
    return(met_NA())
  }
})