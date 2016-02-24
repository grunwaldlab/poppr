library("shiny")
library("poppr")
source("../utils.R")
#------------------------------------------------------------------------------#
# Here, we query the user's R session to find all of the genind and genclone 
# objects
#------------------------------------------------------------------------------#
globals <- get_globals(c("genind", "genclone"))

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  output$selectUI <- renderUI({
    selectInput("dataset", 
                "choose dataset",
                choices = c(globals,
                            "Example: Pinf",
                            "Example: partial_clone",
                            "Example: Aeut",
                            "Example: nancycats",
                            "Example: microbov",
                            "Example: H3N2"),
                selected = "Example: partial_clone"
    )
  })

  #-------------------------------------
  # If the data set is an example, load 
  # and return it, otherwise, get the
  # data set from the user's environment
  #-------------------------------------
  in_dataset <- reactive({
    if (!is.null(input$dataset) && !grepl("<choose>", input$dataset)){
      if(grepl("Example: ", input$dataset)){
        env <- new.env()
        if (input$dataset == "Example: microbov"){ 
          data("microbov", package="adegenet", envir = env) 
        }
        else if (input$dataset == "Example: nancycats"){ 
          data("nancycats", package="adegenet", envir = env) 
        }
        else if (input$dataset == "Example: H3N2"){ 
          data("H3N2", package="adegenet", envir = env) 
        }
        else if (input$dataset == "Example: partial_clone"){ 
          data("partial_clone", package="poppr", envir = env) 
        }
        else if (input$dataset == "Example: Aeut"){ 
          data("Aeut", package="poppr", envir = env) 
        }
        else if (input$dataset == "Example: Pinf"){ 
          data("Pinf", package="poppr", envir = env) 
        }
        exam <- substr(input$dataset, start = 10, stop = nchar(input$dataset))
        dat <- get(exam, envir = env)
      } else {
        dat <- get(input$dataset, envir = .GlobalEnv)
      }
    } else {
      dat <- new("genind")
    }
    if (input$genclone) dat <- as.genclone(dat)
    return(dat)
  })
  
  
  the_method <- reactive({
    input$method
  })
  
  by_pop <- reactive({
    input$bypop
  })
  
  sum2one <- reactive({
    input$sum2one
  })
  
  divisor <- reactive({
    input$divisor
  })
  
  multiplier <- reactive({
    input$mult
  })
  
  eps <- reactive({
    if (!input$epschk) return(NULL)
    if (input$eps > .Machine$double.eps ^0.5){
      return(input$eps)
    } else {
      return(NULL)
    }
  })
  
  pop_arg <- reactive({
    if (input$pop == ""){
      return(NULL)
    }
  })
  
  psex_vector <- reactive({
    input$submit
    poppr::psex(gid = in_dataset(),
                pop = pop_arg(),
                by_pop = by_pop(),
                method = the_method(),
                sum_to_one = sum2one(),
                d = divisor(),
                m = multiplier(),
                e = eps())
  })
  
  output$data_out <- shiny::renderPrint({
    show(in_dataset())
    cat(pop_arg(), by_pop(), the_method(), sep = "\n")
  })
  
  output$psex_plot <- renderPlot({
    input$submit
    if (is.null(psex_vector())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "Please select data and click\nthe 'Go!' button.", 
           cex = 1.6, col = "white")
    } else {
      plot(psex_vector(), log = "y", col = ifelse(psex_vector() > 0.05, "red", "blue"))
      abline(h = 0.05, lty = 2)
    }
  })
  
})
