# Creates a scatterplot comparing MetaboAge and age
output$scatter_metaboage <- renderPlotly({
  if(required()){
    tryCatch({
      req(phenotypes())
      if(length(which(colnames(phenotypes())=="age"))==0){
        plotly_NA_message(main="Age is not available")
      }else{
        x<-data.frame(phenotypes()[,"age"])
        rownames(x)<-rownames(phenotypes())
        scatterplot_predictions(x, MetaboAge(), 
                                title="Chronological Age vs MetaboAge",
                                xname="Chronological age", 
                                yname="MetaboAge")
        }
      }, error = function(err) {
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    })
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})

#Creates an histogram of the mortality score divided for the different ages
output$hist_mort <- renderPlotly({
  if(required()){
      if(dim(phenotypes())[1]!=0){
        if(length(which(colnames(phenotypes())=="age"))==0){
          plotly_NA_message(main="age is not available")
        }else{
          hist_plots_mortality(mort_score(),phenotypes())
        }
      }else{
      return(plotly_NA_message(main="Phenotypes not available,\nplease check your uploaded files."))
    }
  }else{
    return(plotly_NA_message(main="Metabolites not available,\nplease check your uploaded files."))
  }
})


output$pred_choices = renderUI({
  selectInput("pred_choices2",
              "Possible phenotypes:",
              as.character(colnames(predictors()[-1])))
})

output$pheno_choices = renderUI({
    selectInput("pheno_chosen2",
               "Possible phenotypes:",
                as.character(bin_pheno_names()))
})

  
# Outputs the boxplot of the predictor divided by the binary phenotype
output$ttest_score <- renderPlotly({
  tryCatch({
    req(phenotypes())
    pred<-as.character(input$pred_choices2)
    pheno<-as.character(input$pheno_chosen2)
    dat<-data.frame(predictor=predictors()[,pred],pheno=phenotypes()[,pheno])
    colnames(dat)<-c("predictor","pheno")
    ttest_scores(dat = dat, pred=pred,pheno=pheno)
  }, error = function(err) {
    return(plotly_NA_message(main="Image not available,\nplease check your uploaded files."))
  })
})


output$pred_km = renderUI({
  selectInput("pred_km2",
              "Possible predictors:",
              as.character(colnames(predictors()[-1])))
})



# Outputs the boxplot of the predictor divided by the binary phenotype
output$km_score <- renderPlotly({
  if(!is.null(phenotypes())){
    if(length(which(colnames(phenotypes())=="age"))==0){
      plotly_NA_message(main="Age is not available")
    }else if(length(which(colnames(phenotypes())=="Event"))==0 & length(which(colnames(phenotypes())=="EventAge"))==0 ){
      plotly_NA_message(main="Event is not available")
    }else{
      suppressWarnings(kapmeier_scores(predictors=predictors(), pheno=phenotypes(), score=input$pred_km2, Eventname=input$Eventname))
    }
  } else {
    return(plotly_NA_message(main="Image not available,\nplease check your uploaded files."))
  }
})


output$metWASvar = renderUI({
  selectInput("metWASvar2",
              "Test variable:",
              as.character(c(colnames(phenotypes()),
                             colnames(predictors()[-1])
              ))
  )
})

output$metWAScovar = renderUI({
  checkboxGroupInput("metWAScovar2",
                     "Covariates:",
                     as.character(c(colnames(phenotypes()),
                                    colnames(predictors()[-1]))))
})

pred_pheno<-reactive({
  # pheno<-data.frame(ID=rownames(phenotypes()),phenotypes())
  # pred_pheno<-merge(predictors(),pheno, by.x= 'ID',by.y="ID", all=TRUE)
  pred_pheno<-data.frame(predictors(),phenotypes())
  rownames(pred_pheno)<-predictors()$ID
  return(pred_pheno[,-1])
})

# Outputs the boxplot of the predictor divided by the binary phenotype
output$manhplot <- renderPlotly({
  if(!is.null(phenotypes())){
    scaled_met<-data.frame(scale(metabo_measures()))
    rownames(scaled_met)<-rownames(metabo_measures())
    if(!is.null(input$metWAScovar2)){
      phenot<-data.frame(cbind(pred_pheno()[,input$metWASvar2], pred_pheno()[,input$metWAScovar2]))
      colnames(phenot)<-c(input$metWASvar2,input$metWAScovar2)
      rownames(phenot)<-rownames(pred_pheno())
    }else{
      phenot<-cbind(pred_pheno()[,input$metWASvar2])
      colnames(phenot)<-c(input$metWASvar2)
      rownames(phenot)<-rownames(pred_pheno())
    }
    
    if(!is.null(input$metWAScovar2)){
      metWAS<-MetaboWAS(met=scaled_met, pheno=phenot, 
                        test_variable=input$metWASvar2, 
                        covariates=input$metWAScovar2, img=T)
    }else{
      metWAS<-MetaboWAS(met=scaled_met, pheno=phenot, 
                        test_variable=input$metWASvar2, covariates=NULL, img=T)
    }
    
    return(metWAS$manhplot)
  } else {
    return(plotly_NA_message(main="Image not available,\nplease check your uploaded files."))
  }
})


# Outputs the boxplot of the predictor divided by the binary phenotype
res_metWAS <- reactive({
  if(!is.null(phenotypes())){
    scaled_met<-data.frame(scale(metabo_measures()))
    rownames(scaled_met)<-rownames(metabo_measures())
    if(!is.null(input$metWAScovar2)){
      phenot<-data.frame(cbind(pred_pheno()[,input$metWASvar2], pred_pheno()[,input$metWAScovar2]))
      colnames(phenot)<-c(input$metWASvar2,input$metWAScovar2)
      rownames(phenot)<-rownames(pred_pheno())
    }else{
      phenot<-cbind(pred_pheno()[,input$metWASvar2])
      colnames(phenot)<-c(input$metWASvar2)
      rownames(phenot)<-rownames(pred_pheno())
    }
    
    if(!is.null(input$metWAScovar2)){
      metWAS<-MetaboWAS(met=scaled_met, pheno=phenot, 
                        test_variable=input$metWASvar2, 
                        covariates=input$metWAScovar2, img=T)
    }else{
      metWAS<-MetaboWAS(met=scaled_met, pheno=phenot, 
                        test_variable=input$metWASvar2, covariates=NULL, img=T)
    }
    return(metWAS$res[,1:5])
  }
})


## Render MetaboAge to ui
output[["res_metWAS_table"]]<- rendertable(res_metWAS())


