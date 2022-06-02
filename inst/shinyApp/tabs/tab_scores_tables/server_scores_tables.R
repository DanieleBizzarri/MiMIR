# Show results
## Render mortality score to ui
output[["mortScore_table"]]<- rendertable(mort_score())
# output[["mortScore_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(mort_score(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(c(
#       "No data available"
#     )), rownames = FALSE, colnames = ""))
#   })
# })

## Render MetaboAge to ui
output[["MetaboAge_table"]]<- rendertable(MetaboAge())
# output[["MetaboAge_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(MetaboAge(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(c(
#       "No data available"
#     )), rownames = FALSE, colnames = ""))
#   })
# })

## Render T2Dscore to ui
output[["T2Dscore_table"]]<- rendertable(T2D_score_AholaOlli())
# output[["T2Dscore_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(T2D_score_AholaOlli(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(
#       T2D_score_AholaOlli()
#     ), rownames = FALSE, colnames = ""))
#   })
# })

## Render CVDscore to ui
output[["CVDscore_table"]]<- rendertable(CVD_score())
# output[["CVDscore_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(CVD_score(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(
#       CVD_score()
#     ), rownames = FALSE, colnames = ""))
#   })
# })

## Render COVID to ui
output[["COVIDscore_table"]]<- rendertable(COVID_score())
# output[["COVIDscore_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(COVID_score(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(
#       COVID_score()
#     ), rownames = FALSE, colnames = ""))
#   })
# })

## Render MetaboAge to ui
output[["predictors_table"]]<- rendertable(predictors())
# output[["predictors_table"]] <- DT::renderDataTable({
#   tryCatch({
#     DT::datatable(predictors(), options = list(pageLength = 5, scrollX = TRUE))
#   }, error = function(err) {
#     return(DT::datatable(data.frame(
#       CVD_score()
#     ), rownames = FALSE, colnames = ""))
#   })
# })