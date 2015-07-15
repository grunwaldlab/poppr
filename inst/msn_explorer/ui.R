
library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel(paste("Minimum spanning networks in", "poppr")),

  sidebarLayout(
    sidebarPanel(
      tags$h5("status"),
      conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
        tagAppendChild(tags$div(class="progress"),
        tagAppendChild(tags$div(class="progress-bar progress-bar-success", 
                       role="progressbar", `aria-valuenow`="100", 
                       `aria-valuemin`="0", `aria-valuemax`="100", 
                       style="width: 100%"), tags$strong("ready")))
      ),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
        tagAppendChild(tags$div(class="progress"),
        tagAppendChild(tags$div(class="progress-bar progress-bar-striped active", 
                       role="progressbar", `aria-valuenow`="100", 
                       `aria-valuemin`="0", `aria-valuemax`="100", 
                       style="width: 100%"), tags$strong("loading")))
      ),
      tagAppendChildren(
        tags$div(style="display:inline-block"),
        list(
          actionButton("submit", "Go!", icon("check-circle")),
          actionButton("update-data", "reData", icon("refresh")),
          actionButton("update-graph", "reGraph", icon("refresh"))
        )
      ),
      h3("Data Parameters"),
      uiOutput("selectUI"),
      uiOutput("selectPops"),
      
      checkboxInput("genclone", "convert to genclone?", TRUE),
      selectInput("distance", 
                  "choose distance calculation", 
                  choices = c("Dissimilarity",
                              "Bruvo",
                              "Nei",
                              "Rogers",
                              "Edwards",
                              "Provesti",
                              "Reynolds",
                              "Custom")
      ),
      conditionalPanel("input.distance == 'Custom'",
        uiOutput("customDist")
      ),
      conditionalPanel("input.distance == 'Bruvo'",
        selectInput("bruvo_model",
                    "Select a model for missing data",
                    choices = c("Genome Addition",
                                "Genome Loss",
                                "Infinite",
                                "Average Addition/Loss"),
                    selected = "Average Addition/Loss"),
        textInput("replen", "SSR repeat lengths\n(comma separated or a valid R expression)", "1, 2, 3")
      ),
      conditionalPanel("input.distance != 'Bruvo'",
        uiOutput("distargsUI")
      ),
      checkboxInput("reticulate", "include reticulations?", TRUE), 
      h3("Graphical Parameters"),
      checkboxInput("pop.leg", "population legend", TRUE), 
      checkboxInput("scale.leg", "scale bar", TRUE), 
      sliderInput("greyslide",
                  "Grey Scale",
                  min = 0,
                  max = 25,
                  value = 3,
                  step = 1
      ),
      numericInput("nodebase",
                   "Node Size Scale (log(size, value))",
                   "1.15", 
                   min = 1.0001,
                   step = 0.0001),
      numericInput("seed", 
                   "Random Seed",
                   "69"
      ),
      radioButtons("ind_or_mlg", "Labels", 
                   choices = c("sample names", "MLGs"),
                   selected = "sample names", inline = TRUE
      ),
      textInput("inds", NULL, "ALL"),
      checkboxInput("mlgs", "show MLG", FALSE),
      radioButtons("pal", "Indicate a color palette to be used",
                   choices=c("rainbow", 
                            "cm.colors", 
                            "topo.colors", 
                            "terrain.colors", 
                            "gray.colors",
                            "funky",
                            "spectral",
                            "seasun",
                            "azur",
                            "wasp",
                            "custom"), inline = TRUE
      ),
      conditionalPanel("input.pal == 'custom'",
        textInput("custom_pal", "custom palette/function", "function(x) 'purple'")
      ),
      numericInput("cutoff",
                   "Distance Cutoff",
                   NULL,
                   step = 0.001
      ),
      checkboxInput("beforecut", "Keep graph position", TRUE)
    ),


    mainPanel(
      tabsetPanel(
          tabPanel("Plot", plotOutput("plot", height = '600px')),
          tabPanel("Data", verbatimTextOutput("summary")),
          tabPanel("Command", verbatimTextOutput("cmd")),
          tabPanel("Save Plot",
             radioButtons("pdf_png", label = "Choose output filetype",
                          choices = c("pdf", "png"),
                          selected = "pdf",
                          inline = TRUE),
             conditionalPanel("input.pdf_png == 'pdf'",
               numericInput("pdf_plot_width", "Width (in)",
                           value = 7,
                           step = 0.1,
                           min = 1, 
                           max = 20),
               numericInput("pdf_plot_height", "Height (in)",
                           value = 7,
                           step = 0.1,
                           min = 1, 
                           max = 20),
               downloadButton("save_pdf", "Save PDF", class = "btn-info")
            ),
            conditionalPanel("input.pdf_png == 'png'",
              numericInput("png_plot_width", "Width (px)",
                          value = 400,
                          min = 1, 
                          max = 5000),
              numericInput("png_plot_height", "Height (px)",
                          value = 400,
                          min = 1, 
                          max = 5000),
              numericInput("png_res", "Resolution (dpi)",
                          value = 300,
                          min = 72,
                          max = 2000,
                          step = 1),
              downloadButton("save_png", "Save PNG", class = "btn-info")
            )
          ),
          tabPanel("Session Information",
                   verbatimTextOutput("infoRmation")
          )
    )
  )
)
))
