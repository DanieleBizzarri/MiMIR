#Load all the ui of the tabs
source("tabs/tab_upload/ui_upload.R", local = TRUE)
source("tabs/tab_metabolites/ui_metabolites.R", local = TRUE)
source("tabs/tab_phenotypes/ui_phenotypes.R", local = TRUE)
source("tabs/tab_qc/ui_qc.R", local = TRUE)
source("tabs/tab_surro_tables/ui_surro_tables.R", local = TRUE)
source("tabs/tab_scores_tables/ui_scores_tables.R", local = TRUE)
source("tabs/tab_distributions/ui_distributions.R", local = TRUE)
source("tabs/tab_distributions_scores/ui_distributions_scores.R", local = TRUE)
source("tabs/tab_accuracy/ui_accuracy.R", local = TRUE)
source("tabs/tab_accuracy_scores/ui_accuracy_scores.R", local = TRUE)
source("tabs/tab_calibration/ui_calibration.R", local = TRUE)
source("tabs/tab_download/ui_download.R", local = TRUE)
source("tabs/tab_about/ui_about.R", local = TRUE)

#UI for the main menu
ui <- dashboardPage(
  skin = "blue",
  title = "MiMIR",
  
  dashboardHeader(
    # tags$li(class = "dropdown",
    #         tags$style(".main-header {max-height: 20px}"),
    #         tags$style(".main-header .logo {height: 20px;}"),
    #         tags$style(".sidebar-toggle {height: 20px; padding-top: 1px !important;}"),
    #         tags$style(".navbar {min-height:20px !important}"),
            
    #title = span(tagList(icon("laptop-medical"), "MiMIR")),
    title = span(tagList(img(src ="MiMIR_logo_colored_30.jpeg"), "MiMIR")),
                 
    #titlePanel(title=div(img(src="picture.jpg"), "My Title"))
    titleWidth = 350,
    
    
    dropdownMenuOutput("messageMenu"),
    # text with the current tab
    #tags$li(class = "dropdown", tags$a(textOutput("current_page"))),
    #Links to our sites
    tags$li(a(href = 'https://www.lumc.nl/',
              tags$img(src = 'lumc-logo.jpg', height="60", width="60",
                       title = "LUMC"),
              target="_blank"),
            class = "dropdown"),
    tags$li(a(href = 'http://www.molepi.nl/en/home',
              tags$img(src = 'LUMC_MOLEPI_330.png', height="35", width="100",
                  title = "MOLEPI"),
            target="_blank"),
            class = "dropdown"),
    tags$li(a(href = 'https://www.lcbc.nl/',
              tags$img(src = 'logo_lcbc.png', height="60", width="50",
                  title = "LCBC"),
            target="_blank"),
            class = "dropdown"),
    tags$li(a(href = 'https://www.bbmri.nl/',
              tags$img(src = 'logo_BBMRI.jpg', height="45", width="90",
                       title = "BBMRI"),
              target="_blank"),
            class = "dropdown")
    
  ),
  dashboardSidebar(
    collapsed = FALSE,
    width = 350,
    #Sidebar menu with tabs
    sidebarMenu(sidebarMenuOutput("sidebar_tabs"),
                style = "font-size: 15px;")
  ),
  dashboardBody(
    useShinyjs(),
    includeCSS("css/styles.css", local = TRUE),
    
    tabItems(
      tab_upload,
      tab_metabolites,
      tab_phenotypes,
      tab_qc,
      tab_surro_tables,
      tab_scores_tables,
      tab_distributions,
      tab_distributions_scores,
      tab_accuracy,
      tab_accuracy_scores,
      tab_calibration,
      tab_download,
      tab_about
      )
    )
)

