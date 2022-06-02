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
  MOLEPI_LCBC_header(),
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

