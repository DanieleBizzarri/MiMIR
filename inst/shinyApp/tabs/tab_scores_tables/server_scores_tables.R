# Show results
## Render mortality score to ui
output[["mortScore_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(mort_score(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

## Render MetaboAge to ui
output[["MetaboAge_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(MetaboAge(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})


## Render surrogates to ui
output[["surrogates_table"]] <- DT::renderDataTable({
  tryCatch({
    DT::datatable(surrogates(), options = list(pageLength = 5, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

## Render calibrated surrogates to ui
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