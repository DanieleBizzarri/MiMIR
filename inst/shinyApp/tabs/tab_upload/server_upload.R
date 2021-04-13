#######################################
##Basic explanation of a new analysis##
#######################################
output[["welcome_text"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "Welcome to the webtool for calculating the metabolic predictors from raw Nightingale Health 1H-NMR metabolomics data
        developed in the groups of MOLEPI and LCBC (Leiden):mortality score [1], metabolic age ('metaboAge') [2] and surrogate clinical variables [3].
        ",
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      span(
        "Please refer to our manuscripts when using these metabolic biomarkers in your works:",
        style = "text-align: justify; display: block; width: 90%"
      ),
      span(
        "[1] van den Akker Erik B. et al., ‘Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease’, Circ. Genomic Precis. Med., vol. 13, no. 5, pp. 541–547, Oct. 2020, doi: 10.1161/CIRCGEN.119.002610.",
        tags$br(),
        "[2] J. Deelen et al., ‘A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals’, Nat. Commun., vol. 10, no. 1, pp. 1–8, Aug. 2019, doi: 10.1038/s41467-019-11311-9",
        tags$br(),
        "[3] unpublished",
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
    )
  )
})

output[["get_started"]] <- renderUI({
  tagList(
    div(
      br(),
      span(
        "For the webtool to work it will need a file with the metabolic measures in a specific format.
        In case you want to know the accuracy of our models and have them calibrated on your dataset 
        you will need to provide also a file with the phenotypes.",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      span(
        "Metabolic measurements file (necessary) instructions:",
        tags$br(),
        "- The input data should be raw. Don't apply any preprocessing steps.",
        tags$br(),
        "- The columns should represent the metabolites and the rows the samples.",
        tags$br(),
        "- Please use the exact same column names (not necessarily in the same order) as used in the example dataset. 
        Notice that it also requires a columns with unique sample IDs.",
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
        You can check the found column names below.
        The models only use subsets of the metabolites measured by the 
        Nightingale Health platform (see paper). If these metabolites are missing 
        of if they don't have the correct names, an error will be returned.",
        style = "text-align: justify; display: block; width: 95%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 95%"
    )
  )
})

output[["upload_files"]] <- renderUI({
  tagList(
    div(
      span(
        tags$br(),
        "Please, insert the path to your dataset on your local device.",
        tags$br(),
        "Don't worry, your data will be transmitted to our server over an encrypted connection and none of it will be saved nor visualized by us!",
        style = "text-align: justify; display: block; width: 90%"
      ),
      br(),
      style = "background-color: #f5f5f5; border: 1px solid #e3e3e3; width: 90%"
    )
  )
})

###################
##Example dataset##
###################
#load example
# example_met <- reactive({
#   file1 <- "exampleData/syntetic_metabolic_dataset.csv"
#   return(read.table(file=file1, header = FALSE, sep=","))
# })
# 
# example_pheno <- reactive({
#   file2 <- "exampleData/synthetic_phenotypic_dataset.csv"
#   return(read.table(file=file2, header = FALSE, sep=","))
# })


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

###################
## Upload files ##
###################
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
    metabo_measures
    })

#Create data objects
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
  phenotypes
})

bin_phenotypes <- reactive({
  phenotypes<-BMI_LDL_eGFR(phenotypes = phenotypes(),metabo_measures = metabo_measures())
  bin_phenotypes<-binarize_all_pheno(phenotypes())
})

pheno_available <- reactive({
  phen_names<-pheno_names[which(pheno_names %in% colnames(phenotypes()))]
  colnames(phenotypes()[,phen_names])[which(colSums(is.na(phenotypes()[,phen_names]))<nrow(phenotypes()[,phen_names]))]
})

bin_pheno_available <- reactive({
  names(bin_phenotypes())[which(colSums(is.na(bin_phenotypes()))<nrow(bin_phenotypes()))]
})

required<-reactive({
  length(which(MET57 %in% colnames(metabo_measures())))==57
})

## Tables with variables found or not in the files
output$required_met <- renderText({
  req(required())
  "All required metabolites were found!"
})
#Metabolites
output[["found_met"]] <- DT::renderDataTable({
  tryCatch({
    req(input$file_samples$datapath)
    found_met<-data.frame(met=MET57, presence=(MET57 %in% colnames(metabo_measures())))
    found_met<-found_met[order(found_met$presence,decreasing = F),]
    
    DT::datatable(found_met,rownames = F, options = list(pageLength = 10, scrollX = TRUE)) %>% 
                    formatStyle("presence",target = 'row',backgroundColor = styleEqual(c(T, F), c("#B3DE69", "#FB8072"))
                  )
  }, error = function(err) {
    return(DT::datatable(data.frame(c(
      "No data available"
    )), rownames = FALSE, colnames = ""))
  })
})

#Phenotypes
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

# Show analyses files
## Render metabo_measures data to ui
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

## Render count data to ui
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




