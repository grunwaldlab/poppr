library("shiny")
library("poppr")
#------------------------------------------------------------------------------#
# The functions below are those that the server utilizes to parse the data. 
# Most of these allow the server to parse user input. 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Get dist, for example is simply a switch so that the user can select distances
# from a pull-down menu. This is used if the user does not select "Custom".
#------------------------------------------------------------------------------#
get_dist <- function(indist){
  indist <- switch(indist,
      Dissimilarity = "diss.dist",
      Bruvo         = "bruvo.dist",
      Nei           = "nei.dist",
      Rogers        = "rogers.dist",
      Edwards       = "edwards.dist",
      Provesti      = "provesti.dist",
      Reynolds      = "reynolds.dist"
  )
  return(indist)
}
#------------------------------------------------------------------------------#
# If the user does select custom, this function will make sure that it is 
# encapsulated in parentheses. This makes sure that the entire expression is
# evaluated.
#------------------------------------------------------------------------------#
parse_distfun <- function(x){
  if (grepl("function", x)){
    x <- paste0("(", x, ")")
  }
  return(x)
}
#------------------------------------------------------------------------------#
# This is a function that will print vectors of numerics or characters in a way
# that can be run directly from R. For example, if the user checks two boxes
# labeling the populations "pop1" and "pop2", that vector of populations gets 
# passed onto the popsub function. By default, it will print like so:
#
# [1] "pop1" "pop2"
#
# You cannot copy and paste this into R because it will throw an error, thus, 
# this function will print the vector like so:
#
# c("pop1", "pop2")
#
# This is then usable in the Command tab for the popsub function.
#------------------------------------------------------------------------------#
make_dput <- function(x){
  return(capture.output(dput(x)))
}

#------------------------------------------------------------------------------#
# A function to search the user's global environment and grab all of the useable
# objects. In this case, it's genind and genclone objects, but the objclass
# argument allows this to be extensible to any class. This is immensely useful
# so that the user does not have to save their objects as rda files, nor do they
# have to save them as text files for input.
#------------------------------------------------------------------------------#
get_globals <- function(objclass = c("genind", "genclone")){
  # Grab all object names in users R session.
  myobjs <- ls(envir = .GlobalEnv) 
  if (length(myobjs) == 0){
    return(myobjs)
  }
  # Go through each name and test if it is any of the classes in objclass.
  gens <- vapply(myobjs, FUN = is_usable, FUN.VALUE = logical(1), objclass)
  myobjs[gens]
}

#------------------------------------------------------------------------------#
# This function tests a single object name to see if it is of class objclass.
# The function is used in get_globals
#------------------------------------------------------------------------------#
is_usable <- function(object, objclass = c("genind", "genclone")){
  any(objclass %in% class(get(object, .GlobalEnv)))
}
#------------------------------------------------------------------------------#
# Here, we query the user's R session to find all of the genind and genclone 
# objects
#------------------------------------------------------------------------------#
globals <- get_globals(c("genind", "genclone"))

shinyServer(function(input, output, session) {
  
  
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
  
  output$selectPops <- renderUI({
    input$dataset
    checkboxGroupInput("sublist",
                "choose populations",
                choices = popNames(in_dataset()),
                inline = TRUE,
                selected = popNames(in_dataset()))
  })
  
  dataset <- reactive({
    input$`update-data`
    input$dataset
    isolate({
      popsub(in_dataset(), input$sublist, drop = FALSE)
    })
  })

  dataname <- reactive({
    if (!grepl("<choose>", input$dataset)){
      if(grepl("Example: ", input$dataset)){
        dat <- substr(input$dataset, start = 10, stop = nchar(input$dataset))
      } else {
        dat <- input$dataset
      }
    } else {
      dat <- "no data"
    }
    return(dat)
  })
  
  distfun <- reactive({ 
    if (input$distance == "Custom"){
      parse_distfun(input$custom_distance)
    } else {
      get_dist(input$distance) 
    }
  })

  output$customDist <- renderUI({
    textInput("custom_distance", label = "Custom Distance Function", "function(x) dist(tab(x))")
  })
  
  output$distargsUI <- renderUI({
    the_fun <- eval(parse(text = distfun()))
    the_args <- formals(the_fun)[-1]
    the_args <- paste(names(the_args), unlist(the_args), sep = " = ", 
                      collapse = ", ")
    textInput("distargs", label = "Distance arguments", the_args)
  })

  reticulation <- reactive({
    input$reticulate        
  })
  
  
  distargs <- reactive({
    input$distargs     
  })
  
  addloss <- reactive({
    switch(input$bruvo_model,
           "Genome Addition" = "add = TRUE, loss = FALSE",
           "Genome Loss" = "add = FALSE, loss = TRUE",
           "Infinite" = "add = FALSE, loss = FALSE",
           "Average Addition/Loss" = "add = TRUE, loss = TRUE")
  })

  replen <- reactive({
    if (!grepl("\\(", input$replen)){
      paste0("replen = c(", input$replen, ")")      
    } else {
      paste0("replen = ", input$replen)
    }
  })

  minspan <- reactive({
    input$dataset
    input$`update-data`
    input$reticulate
    input$submit
    isolate({
      indist <- distfun()
      ret    <- reticulation()
      args   <- distargs()
      if (input$distance == "Bruvo"){
        args <- paste(replen(), addloss(), sep = ", ")
        fun <- paste0("bruvo.msn(dataset(), ", args, ", showplot = FALSE, include.ties = ret)")
        out <- eval(parse(text = fun))
      } else {
        if (indist != "diss.dist"){
          dat <- missingno(dataset(), "mean")
        } else {
          dat <- dataset()
        }
        if (length(args) == 1 && args == ""){
          fun <- paste0(indist, "(dat)")
        } else {
          fun <- paste0(indist, "(dat, ", args, ")")
        }
        dist <- eval(parse(text = fun))
        out <- poppr.msn(dataset(), dist, showplot = FALSE, include.ties = ret)
      }
      return(out)
    })
  })

  slide <- reactive({
    # isolate({
      return(input$greyslide)
    # })
  })
  
  seed <- reactive({
    # isolate({
      return(input$seed)      
    # })

  })

  nodebase <- reactive({
    # isolate({
      return(input$nodebase)
    # })
  })

  inds <- reactive({
    # isolate({
      inds <- strsplit(input$inds, "[[:blank:]]*,[[:blank:]]*")[[1]]

      if (input$ind_or_mlg == "sample names" || inds == "ALL" || inds == ""){
        return(inds)
      } else {
        return(as.numeric(inds))
      }
    # })

  })
  
  usrPal <- reactive({
    input$`update-data`
    input$`update-graph`
    isolate({
      if (input$pal == 'custom'){
        eval(parse(text = input$custom_pal))
      } else {
        input$pal
      }
    })
  })

  popLeg <- reactive({
    # isolate({
      input$pop.leg
    # })
  })

  scaleLeg <- reactive({
    # isolate({
      input$scale.leg
    # })
  })

  cutoff <- reactive({
    cutoff <- as.numeric(input$cutoff)
    if (is.na(cutoff)) cutoff <- NULL
    cutoff      
  })
  
  distcmd <- reactive({
    dat      <- dataname()
    distfunk <- distfun()
    args     <- distargs()
    the_pops <- popNames(in_dataset())
    match_pops <- the_pops %in% input$sublist
    half <- ceiling(length(the_pops)/2)
    if (sum(match_pops) < half){
      first_dat <- paste0(dat, "_sub <- popsub(", dat, ", sublist = ", make_dput(input$sublist), ")\n")
    } else {
      first_dat <- paste0(dat, "_sub <- popsub(", dat, ", blacklist = ", make_dput(the_pops[!match_pops]), ")\n")
    }
    closer   <- paste0("showplot = FALSE, include.ties = ", reticulation(), ")")
    has_no_args <- length(args) == 1 && args == ""
    if (distfunk == "bruvo.dist"){
      args <- paste(replen(), addloss(), sep = ", ")
      distfunk <- "min_span_net <- bruvo.msn"
      closer <- paste0(", ", args, ", ", closer)
      return_cmd <- paste0(distfunk, "(", dat, "_sub", closer)
    } else { 
      if (distfunk == "diss.dist"){
        missfunk <- character(0)
        distfunk <- paste0(distfunk, "(", dat, "_sub, ", args, ")\n")        
      } else {
        missfunk <- paste0(dat, "_nomiss <- ", "missingno(", dat, 
                           ", type = 'mean')\n")
        args <- ifelse(has_no_args, "", paste0(", ", args))
        distfunk <- paste0(distfunk, "(", dat, "_nomiss", args, ")\n")        
      }
      msnfunk <- paste0("poppr.msn(", dat, "_sub, ", dat, "_dist, ", closer, "\n")
      return_cmd <- paste0(missfunk, 
                           dat, "_dist <- ", distfunk,
                           "min_span_net <- ", msnfunk)
    }
    return(paste0(first_dat, return_cmd))
  })
  
  cmd <- reactive({
    dat <- dataname()
    pal <- ifelse(input$pal == 'custom', input$custom_pal, input$pal)
    paste0("plot_poppr_msn(", dat, 
           ",\n\t       min_span_net", 
           ",\n\t       inds = ", make_dput(inds()), 
           ",\n\t       mlg = ", input$mlgs,
           ",\n\t       gadj = ", input$greyslide,
           ",\n\t       nodebase = ", input$nodebase,
           ",\n\t       palette = ", pal,
           ",\n\t       cutoff = ", ifelse(is.null(cutoff()), "NULL", cutoff()),
           ",\n\t       quantiles = FALSE",
           ",\n\t       beforecut = ", bcut(), ")")
  })

  bcut <- reactive({
    input$beforecut
  })

  output$summary <- renderPrint({
    dat <- dataset()
    show(dat)
  })

  output$plot <- renderPlot({
    input$pop.leg
    input$scale.leg
    input$beforecut
    input$nodebase
    input$inds
    input$mlgs
    input$`update-graph`
    if(!input$submit) {
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "Please select data and click\nthe 'Go!' button.", 
           cex = 1.6, col = "white")
    } else {
      set.seed(seed())
      plot_poppr_msn(dataset(), 
                     minspan(), 
                     ind = inds(), 
                     gadj = slide(), 
                     mlg = input$mlgs,
                     palette = usrPal(), 
                     cutoff = cutoff(), 
                     quantiles = FALSE, 
                     beforecut = bcut(), 
                     nodebase = nodebase(),
                     pop.leg = popLeg(), 
                     scale.leg = scaleLeg()
                    )      
    }
  })
  
  output$cmd <- renderPrint({
    cat(paste0(distcmd(), "\n"))
    cat(paste0("set.seed(", seed(),")\n"))
    cat(cmd())
  })

  output$save_pdf <- downloadHandler(
    filename = function() paste0('msn-', Sys.Date(), '.pdf'),
    content = function(file) {
      isolate({
        # Generate a png
        pdf(file, width = input$pdf_plot_width, height = input$pdf_plot_height)
        set.seed(seed())
        plot_poppr_msn(dataset(),
                       minspan(),
                       ind = inds(),
                       gadj = slide(),
                       mlg = input$mlgs,
                       palette = usrPal(),
                       cutoff = cutoff(),
                       quantiles = FALSE, 
                       beforecut = bcut(),
                       nodebase = nodebase(),
                       pop.leg = popLeg(),
                       scale.leg = scaleLeg()
                      )
        dev.off()
      })      
    }
  )
  output$save_png <- downloadHandler(
    filename = function() paste0('msn-', Sys.Date(), '.png'),
    content = function(file) {
      isolate({
        # Generate a png
        png(file, width = input$png_plot_width, height = input$png_plot_height)
        set.seed(seed())
        plot_poppr_msn(dataset(),
                       minspan(),
                       ind = inds(),
                       gadj = slide(),
                       mlg = input$mlgs,
                       palette = usrPal(),
                       cutoff = cutoff(),
                       quantiles = FALSE, 
                       beforecut = bcut(),
                       nodebase = nodebase(),
                       pop.leg = popLeg(),
                       scale.leg = scaleLeg()
                      )
        dev.off()
      })      
    }
  )
  output$infoRmation <- renderPrint({
    sessionInfo()
  })
})
