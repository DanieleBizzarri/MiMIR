tab_about <- tabItem(
  tabName = "about",
  align = "center",
  
  h1("MetaboRiSc"),
  uiOutput("description"),
  # h2("Models description"),
  # tags$img(src = "Intro_plot.png", height = 400, width = 700),
  HTML('<hr style="border-color: #0088cc;">'),
  
  h3("Version: 1.0"),
  h2("Developed by the LCBC team (LUMC)"),
  h3("Daniele Bizzarri", tags$br(),
     "Erik van den Akker", tags$br(),
     "Marcel Reinders", tags$br(),
     "Eline Slagboom", tags$br()),
  br(),
  HTML('<hr style="border-color: #0088cc;">'),
  
  h2("Session info"),
  verbatimTextOutput("currentSession")
  
  # HTML('<hr style="border-color: #0088cc;">'),
  # tags$img(src = "LUMC_MOLEPI_330.png", height = 60, width = 200),
  # tags$img(src = "lcbc_logo2.png", height = 60, width = 250),
  # tags$img(src = "logo_BBMRI.jpg", height = 60, width = 150)
  
  # tags$li(a(href = 'http://www.molepi.nl/en/home',
  #           tags$img(src = "LUMC_MOLEPI_330.png", height = 60, width = 200,
  #                    title = "MOLEPI"),
  #           target="_blank"),class = "dropdown"),
  # tags$li(a(href = 'https://www.lcbc.nl/',
  #           tags$img(src = "lcbc_logo2.png", height = 60, width = 250,
  #                    title = "LUMC"),
  #           target="_blank"),class = "dropdown"),
  # tags$li(a(href ='https://www.bbmri.nl/',
  #           tags$img(src = 'logo_BBMRI.jpg', height = 60, width = 200,
  #                    title = "BBMRI-NL"),
  #           target="_blank"),class = "dropdown"),
  
)

