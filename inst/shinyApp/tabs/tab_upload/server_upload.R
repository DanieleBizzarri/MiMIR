#Basic explanation to start new analysis
output[["get_started"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "This webtool was created to easily allow the computation of several
        existing metabolic scores trained using Nightingale Health data.
        This application works only with file, including the Nighitngale NH-metabolomics measurements, in a specific format:",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      span(
        "Metabolic measurements file (necessary) instructions:",
        tags$br(),
        "- Input data should be raw. Don't apply any preprocessing steps.",
        tags$br(),
        "- The columns: metabolic features and rows samples.
        Please include also the column names with the metabolites names and an 'ID' column with unique sample IDs.",
        tags$br(),
        "- The app can translate alternative metabolic names created by Nightingale health.
        The table 'Variables names recognised', below, will allow you to identify unrecognised features.
        If a metabolite is not found please update your file with the name you can find in the example dataset.
        The models only use subsets of the metabolites measured by the 
        Nightingale Health platform (see Supplementary Info). If these metabolites are missing 
        or if they don't have the correct names, an error will be returned.",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      span(
        "It is also recommended to provide a file with the required phenotypes to calculate the accuracy and
        calibrate the models. Phenotypic variables file (optional) instructions:",
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
        You can check the column names below in the table 'Variables names recognised'. 
        If a column updated was not recognised please upload again the file with the correct nomenclature",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 95%"
    )
  )
})


####################
## Uploaded files ##
####################
## Read metabolites data file
upload_met <- reactive({
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
    #avoid case-sensitive alternative names
    colnames(metabo_measures)<-tolower(colnames(metabo_measures))
    #Looking for alternative names
    nam<-find_BBMRI_names(colnames(metabo_measures))
    i<-which(nam$BBMRI_names %in% MiMIR::metabo_names_translator$BBMRI_names)
    metabo_measures<-metabo_measures[,i]
    colnames(metabo_measures)<-nam$BBMRI_names[i]
    
    if(length(which(colnames(metabo_measures)=="faw6_faw3"))==0){
      metabo_measures$faw6_faw3<-metabo_measures$faw6/metabo_measures$faw3
    }
    metabo_measures
    })

#Read phenotypes data file
upload_phen <- reactive({
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
  
  if(!is.null(metabo_measures())){
    phenotypes<-phenotypes[rownames(metabo_measures()),]
  }
  
  #avoid case-sensitive alternative names
  #colnames(phenotypes)<-tolower(colnames(phenotypes))
  
  phenotypes
})

####################
## Synthetic Data ##
####################
syn_met<-eventReactive(input$load_synth_data,{
  metabo_measures<-as.data.frame(MiMIR::synthetic_metabolic_dataset)
})

syn_phen<-eventReactive(input$load_synth_data,{
  phenotypes<-MiMIR::synthetic_phenotypic_dataset
})

# Download the example synthetic dataset
shiny::observeEvent(
  eventExpr = input$downloadData,
  handlerExpr = {
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    session = getDefaultReactiveDomain()
    shinyFiles::shinyFileSave(input = input,id = "downloadData",roots = volumes,session = session,
                              restrictions = system.file(package = "base"))
    fileinfo <- shinyFiles::parseSavePath(roots = volumes, selection = input$downloadData)
    if (nrow(fileinfo) > 0){
      temp <- setwd(tempdir())
      setwd(temp)
      files <- c("example_metabolic_dataset.csv", "example_phenotypic_dataset.csv")
      
      write.csv(MiMIR::synthetic_metabolic_dataset, "example_metabolic_dataset.csv")
      write.csv(MiMIR::synthetic_phenotypic_dataset, "example_phenotypic_dataset.csv")
      zip(zipfile = fileinfo$datapath, files = files)
    }
  }
)

####################################
## Dataset (uploaded or synthetic)##
####################################
metabo_measures<-reactive({
  if(!is.null(upload_met())){
    metabo_measures<-upload_met()
  }else if(!is.null(syn_met())){
    metabo_measures<-syn_met()
  }else{
    metabo_measures<-as.data.frame(MiMIR::synthetic_metabolic_dataset)
  }
})

phenotypes<-reactive({
  if(!is.null(upload_phen())){
    phenotypes<-upload_phen()
  }else if(!is.null(syn_met())){
    phenotypes<-syn_phen()
  }else{
    phenotypes<-data.frame()
  }
})

##############################
## Check available features ##
##############################
# Variable storing the selected phenotypes found in the uploaded file
pheno_available <- reactive({
  phen_names<-pheno_names[which(phenotypes_names$pheno_names %in% colnames(phenotypes()))]
  colnames(phenotypes()[,phen_names])[which(colSums(is.na(phenotypes()[,phen_names]))<nrow(phenotypes()[,phen_names]))]
})

bin_pheno_names<-reactive({
  #colnames(phenotypes())
  colnames(phenotypes())[apply(phenotypes(),2,function(x) {all(na.omit(x) %in% 0:1) })]
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
  if(length(metabo_measures())>0){
    length(which(MiMIR::metabolites_subsets$MET14 %in% colnames(metabo_measures())))==14
  }else{FALSE}
})

## Tables with variables found or not in the files
output$required_met <- renderText({
  if(required()){
    "All required metabolites were found!"
  }else{"Some required metabolites were not found!"}
})

############
## Tables ##
############
#Table with the metabolites names found in the uploaded file
output[["found_met"]] <- DT::renderDataTable({
  tryCatch({
    req(metabo_measures())
    found_met<-data.frame(met=MiMIR::metabolites_subsets$MET62, presence=(MiMIR::metabolites_subsets$MET62 %in% colnames(metabo_measures())))
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
    req(phenotypes())
    found_phen<-data.frame(phen=MiMIR::phenotypes_names$pheno_names, presence=(MiMIR::phenotypes_names$pheno_names %in% colnames(phenotypes())))
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
    req(metabo_measures())
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
    req(phenotypes())
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
    req(phenotypes())
    DT::datatable(bin_phenotypes(), options = list(pageLength = 10, scrollX = TRUE))
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})




