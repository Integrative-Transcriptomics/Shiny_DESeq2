source("./Components/ShinySeq_Packages.R")
source("./Components/ShinySeq_Functions.R")
source("./Components/ShinySeq_UI.R")
source("./Components/ShinySeq_Server.R")

options(shiny.maxRequestSize=200*1024**2) 
shinyApp(ui, server)

