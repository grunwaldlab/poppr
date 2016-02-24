#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Show Psex values"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("Data Input"),
      uiOutput("selectUI"),
      checkboxInput("genclone", "convert to genclone?", TRUE),
      h3("Psex calculations"),
      checkboxInput("bypop", "by population", TRUE),
      textInput("pop", "strata", 
                value = "", 
                placeholder = "~Population/Subpopulation"),
      radioButtons("method", "Method",
                   choices = c("single", "multiple"),
                   selected = "single", inline = TRUE),
      h3("Minor Allele Frequency Correction"),
      checkboxInput("epschk", "Custom correction", FALSE),
      sliderInput("eps", "Correct all minor allele frequencies to:", 1e-5,
                   min = 1e-5, max = 1, step = 0.0001),
      checkboxInput("sum2one", "sum allele frequencies to one?", FALSE),
      radioButtons("divisor", "Divisor",
                   choices = c("sample", "mlg", "rrmlg"),
                   selected = "sample", inline = TRUE),
      sliderInput("mult", "Multiplier", 1, min = 0.001,
                   max = 1, step = 0.05),
      actionButton("submit", "Go!", icon("check-circle"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      shiny::plotOutput("psex_plot")
    )
  )
))
