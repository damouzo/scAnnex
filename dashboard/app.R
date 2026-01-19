# scAnnex Dashboard - Main Application
# Entry point for the Shiny app

# Source global functions and UI/Server
source("global.R")
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)
