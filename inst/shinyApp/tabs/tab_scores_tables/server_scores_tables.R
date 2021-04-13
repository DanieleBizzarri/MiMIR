# Show results
## Render mort_score data to ui
output[["mortScore_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(mort_score(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

## Render MetaboAge data to ui
output[["MetaboAge_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(MetaboAge(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})


## Render mort_score data to ui
output[["surrogates_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(surrogates(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

# calibrations_data.frame<-reactive({
#   cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
# })

output[["calibrated_surro_table"]] <- DT::renderDataTable({
  cal<-calib_data_frame(calibrations(), metabo_measures(), bin_pheno_available())
  tryCatch({
    DT::datatable(data.frame(cal), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})