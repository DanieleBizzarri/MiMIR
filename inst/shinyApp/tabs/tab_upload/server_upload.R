#Basic explanation to start new analysis
output[["get_started"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "This webtool was created with the intention of providing a set of interesting scores
        retrieved from the metabolic features created by Nightingale-Health. 
        Thee app works only with a file with the Nighitngale NH-metabolomics measurements  in a specific format.
        It is also recommended to provide a file with the required phenotypes to calculate the accuracy and
        calibrate the surrogate models.",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      span(
        "Metabolic measurements file (necessary) instructions:",
        tags$br(),
        "- Input data should be raw. Don't apply any preprocessing steps.",
        tags$br(),
        "- The columns should represent the metabolites and the rows the samples.",
        tags$br(),
        "- The app can translate alternative metabolic names created by Nightingale health.
        To know if a metabolite wasn't recognised look below at the table 'Variables names recognised'.
        If a metabolite is not found please upload again the file with the name you can find in the example dataset.
        Notice that it also requires a columns with unique sample IDs.
        The models only use subsets of the metabolites measured by the 
        Nightingale Health platform (see paper). If these metabolites are missing 
        of if they don't have the correct names, an error will be returned.",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      span(
        "Phenotypic variables file (optional) instructions:",
        tags$br(),
        "- The input data should contain the raw phenotypes with the columns name as the example dataset.",
        tags$br(),
        "- The file should contain a column with the same IDs as the metabolic measurement file.",
        tags$br(),
        "- It's not necessary that all the phenotypes listed are available!",
        tags$br(),
        style = "text-align: justify; display: flex; block; width: 95%"
      ),
      span(
        tags$br(),
        "Please ensure that your dataset looks similar to this example dataset before submitting.
        You can check the found column names below in the table 'Variables names recognised'. 
        If a column updated was not recognise please upload agein the file with the correct naming",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 95%"
    )
  )
})

###################
##Example dataset##
###################
output$downloadData <- downloadHandler(
  filename = function() {
    paste("metabolic_predictors_example_dataset", "zip", sep=".")
  },
  content = function(fname) {
    temp <- setwd(tempdir())
    setwd(temp)
    files <- c("example_metabolic_dataset.csv", "example_phenotypic_dataset.csv")
    
    write.csv(synthetic_metabolic_dataset, "example_metabolic_dataset.csv")
    write.csv(synthetic_phenotypic_dataset, "example_phenotypic_dataset.csv")
    zip(zipfile = fname, files = files)
  },
  contentType = "application/zip"
  )

####################
## Uploaded files ##
####################
## Read metabolites data file
metabo_measures <- reactive({
  if (is.null(input$file_samples$datapath)) {
    return(NULL)
  }
    metabo_measures <- read.csv(
      input$file_samples$datapath,
      header = TRUE, row.names = 1,
      sep = "\t"
    )
    if (ncol(metabo_measures) <= 1) {
      metabo_measures <- read.csv(
        input$file_samples$datapath,
        header = TRUE, row.names = 1
      )
    }
    
    #Looking for alternative names
    nam<-find_BBMRI_names(colnames(metabo_measures))
    i<-which(nam$BBMRI_names %in% metabo_names_translator$BBMRI_names)
    metabo_measures<-metabo_measures[,i]
    colnames(metabo_measures)<-nam$BBMRI_names[i]
    
    if(length(which(colnames(metabo_measures)=="faw6_faw3"))==0){
      metabo_measures$faw6_faw3<-metabo_measures$faw6/metabo_measures$faw3
    }
    
    metabo_measures
    })

#Read phenotypes data file
phenotypes <- reactive({
  if (is.null(input$file_phenotypes$datapath)) {
    return(NULL)
  }
  phenotypes <- read.csv(
    input$file_phenotypes$datapath,
    header = TRUE, row.names = 1,
    sep = "\t"
  )
  if (ncol(phenotypes) <= 1) {
    phenotypes <- read.csv(
      input$file_phenotypes$datapath,
      header = TRUE,row.names = 1
    )
  }
  
  # if(!is.null(metabo_measures())){
  #   ind<-order(match(rownames(phenotypes()), rownames(metabo_measures())))
  #   predictors<-predictors[ind,]
  # }
  
  if(!is.null(metabo_measures())){
    phenotypes<-phenotypes[rownames(metabo_measures()),]
  }
  
  phenotypes
})

# Variable storing the selected phenotypes found in the uploaded file
pheno_available <- reactive({
  phen_names<-pheno_names[which(pheno_names %in% colnames(phenotypes()))]
  colnames(phenotypes()[,phen_names])[which(colSums(is.na(phenotypes()[,phen_names]))<nrow(phenotypes()[,phen_names]))]
})

bin_pheno_names<-reactive({
  #colnames(phenotypes())
  colnames(phenotypes())[apply(phenotypes(),2,function(x) { all(na.omit(x) %in% 0:1) })]
})

#Create the binarize phenotypes table
bin_phenotypes <- reactive({
  phenotypes<-BMI_LDL_eGFR(phenotypes = phenotypes(),metabo_measures = metabo_measures())
  bin_phenotypes<-binarize_all_pheno(phenotypes())
  
  ind<-order(match(rownames(bin_phenotypes), rownames(metabo_measures())))
  bin_phenotypes<-bin_phenotypes[ind,]
  return(bin_phenotypes)
})

# Variable storing the selected phenotypes found in the binarized phenotypes table
bin_pheno_available <- reactive({
  names(bin_phenotypes())[which(colSums(is.na(bin_phenotypes()))<nrow(bin_phenotypes()))]
})

# Variable TRUE/FALSE if all the metabolites names were found
required<-reactive({
  length(which(MET57 %in% colnames(metabo_measures())))==57
})

## Tables with variables found or not in the files
output$required_met <- renderText({
  req(required())
  "All required metabolites were found!"
})

#Table with the metabolites names found in the uploaded file
output[["found_met"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_samples$datapath)
    found_met<-data.frame(met=MET62, presence=(MET62 %in% colnames(metabo_measures())))
    found_met<-found_met[order(found_met$presence, decreasing = F),]
    if(found_met[which(found_met$met=="faw6_faw3"),"presence"]==F){
      found_met[which(found_met$met=="faw6_faw3"),"presence"]<-NA
      
    }
    
    DT::datatable(found_met, rownames = F, options = list(pageLength = 10, scrollX = TRUE)) %>% 
                    formatStyle("presence",target = 'row',backgroundColor = styleEqual(c(T, F, NA), c("#B3DE69", "#FB8072",  "#FDB462"))
                  )
    
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

#Table with phenotypes names found in the uploaded file
output[["found_phen"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_phenotypes$datapath)
    found_phen<-data.frame(phen=pheno_names, presence=(pheno_names %in% colnames(phenotypes())))
    found_phen<-found_phen[order(found_phen$presence,decreasing = F),]
    
    DT::datatable(found_phen,rownames = F, options = list(pageLength = 10, scrollX = TRUE)) %>% 
      formatStyle("presence",target = 'row',backgroundColor = styleEqual(c(T, F), c("#B3DE69", "#FDB462"))
      )
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

# Show uploaded files
# Render metabo_measures data to ui
output[["metabo_table"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_samples$datapath)
    DT::datatable(metabo_measures(), options = list(pageLength = 10, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

# Render pheno_table data to ui
output[["pheno_table"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_phenotypes$datapath)
    DT::datatable(phenotypes(), options = list(pageLength = 10, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

# Render bin_pheno_table data to ui
output[["bin_pheno_table"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_phenotypes$datapath)
    DT::datatable(bin_phenotypes(), options = list(pageLength = 10, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})




