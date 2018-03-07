#### ui.R

setwd("E:/Users/cwang/DMGE/")
# setwd("/home/cwang/project/dmge")

ui <- shinyUI(fluidPage(
  # Application title
  headerPanel("CSRE browser V0.3"),
  fluidRow(
    column(4, verbatimTextOutput("region", TRUE)),
    column(4, textInput("Gene_OR_Region", NULL, "",
                        placeholder = "enter a gene symbol or region"),
           helpText("A gene symbol (e.g. PRAME) or region (e.g. chr22:22888118-22902174)")),
    column(4, actionButton("go", "go"))
  ),
  fluidRow(
    column(12, align = "center",
           span("move"),
           span(title = "move 95% to the left",
                actionButton("left_95", "<<<")),

           span(title = "move 47.5% to the left",
                actionButton("left_47.5", "<<")),

           span(title = "move 10% to the left",
                actionButton("left_10", "<")),

           span(title = "move 10% to the right",
                actionButton("right_10", ">")),

           span(title = "move 47.5% to the right",
                actionButton("right_47.5", ">>")),

           span(title = "move 95% to the right",
                actionButton("right_95", ">>>")),

           span("zoom in"),
           actionButton("in_1.5", "1.5x"),
           actionButton("in_3", "3x"),
           actionButton("in_10", "10x"),

           span("zoom out"),
           actionButton("out_1.5", "1.5x"),
           actionButton("out_3", "3x"),
           actionButton("out_10", "10x"))),
  # Show a plot of the CSRE
  uiOutput("plotUI"),
  # div(plotOutput("plotCSRE", height = 400,
  #            hover = hoverOpts("plot_hover", 100, "debounce"),
  #            dblclick = "plot_dblclick",
  #            brush = brushOpts(id = "plot_brush",
  #                              direction = "x",
  #                              resetOnNew = TRUE))),
  plotOutput("highlight", height = 100),
  fluidRow(
    column(4, verbatimTextOutput("subRegion")),
    column(4, verbatimTextOutput("pos"))),
  # Show the table of detail information
  DT::dataTableOutput("ft"),
  downloadButton('downloadData', 'Download this table')
))
