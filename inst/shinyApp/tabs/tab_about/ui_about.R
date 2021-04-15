tab_about <- tabItem(
  tabName = "about",
  align = "center",
  
  h1("MetaboRiSc"),
  #Description of our papers
  uiOutput("description"),
  HTML('<hr style="border-color: #0088cc;">'),
  
  #Description of us
  h3("Version: 1.0"),
  h2("Developed by the LCBC team (LUMC)"),
  h3("Daniele Bizzarri", tags$br(),
     "Erik van den Akker", tags$br(),
     "Marcel Reinders", tags$br(),
     "Eline Slagboom", tags$br()),
  br(),
  HTML('<hr style="border-color: #0088cc;">'),
  
  #Session info
  h2("Session info"),
  verbatimTextOutput("currentSession")
  
)

