suppressMessages(library(shiny))
suppressMessages(library(shinyBS))
suppressMessages(library(shinyLP))
suppressMessages(library(shinythemes))
suppressMessages(library(shinyjs))
suppressMessages(library(shinycssloaders))
suppressMessages(library(DT))
suppressMessages(library(data.table))
suppressMessages(library(plotly))


file_path="final"

#gene_list <- readRDS(file = paste(file_path,"uniprot_genelist.rds",sep="/"))

# Define UI for random distribution application 
shinyUI( fluidPage(

  div(style="padding: 0px 0px; width: '100%'", titlePanel( title="", windowTitle="VarViz" )),
  tagList(
    useShinyjs(),
  navbarPage(title="VarViz",
             inverse = F, 
             theme = shinytheme("united"),
             tabPanel("Home", icon = icon("home"),
                      
                      jumbotron("Hello Geneticist!", "Call attention to important application features or provide guidance", buttonLabel = "Click Me"),
                      fluidRow(
                        column(6, panel_div(class_type = "primary", panel_title = "Directions", content = "How to use the app")),
                        column(6, panel_div("success", "Application Maintainers",
                                            HTML("Email Me: <a href='mailto:ametpally@raleighcharterhs.org?Subject=Gene VarINFO plot%20Help' target='_top'>Agasthya Metpally</a>")))
                      ),  # end of fluidRow
                      fluidRow(
                        column(6, panel_div("info", "App Status", "Include text with status, version and updates")),
                        column(6, panel_div("danger", "Security and License", "Copyright 2024")),
                        
                        #### FAVICON TAGS SECTION ####
                        tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
                        
                        bsModal("modalExample", "Instructional Video", "tabBut", size = "large" ,
                                p("Additional text and widgets can be added in these modal boxes. Video plays in chrome browser"),
                                iframe(width = "560", height = "315", url_link = "https://www.youtube.com/embed/0fKg7e37bQE")
                        )
                        
                      )),
             tabPanel("Protein View",
                      shinyjs::useShinyjs(),
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          uiOutput("selects1"),
                          checkboxGroupInput("plotselection", "Output Plot:",
                                             choices = c("AlphaFold" = "alphafold", 
                                                         "Clinvar" = "clinvar",
                                                         "Density" = "density",
                                                         "gnomAD Freq" = "freq", 
                                                         "gnomAD Depth" = "depth"),
                                             selected = c("alphafold", "clinvar", "density", "freq", "depth")),
                          selectizeInput("gene_name", "Choose a gene", choices = NULL),
                          #selectInput("gene_name", "Choose a Gene:", c(Choose='', gene_list), selectize=TRUE),
                          fileInput("consurf_file", accept = c("text/csv","text/comma-separated-values,text/plain",".txt"), label = h5("Please provide Consurf scores in tab separated format:")),
                          fileInput("user_file", accept = c("text/csv","text/comma-separated-values,text/plain",".txt"), label = h5("Please provide protein locations:")),
                          
                          textAreaInput("variants", "Please provide comma separated variants in pMut format:", "p.Val109Asp", height = "100px"),
                          radioButtons("label","Variant Label:",inline=TRUE, choices = list("Yes","No"), selected = "Yes"),
                          radioButtons("clinvar_filter","Clinvar Variants Filter :", inline=TRUE, choices = list("Any Star","1 Star or more","2 Stars or more"), selected = "1 Star or more"),
                          htmlOutput("maxAC_F"),
                          br(),
                          a(id = "toggleAdvanced", "Show/Hide this section", href = "#"),br(), 
                          shinyjs::hidden(
                            div(id = "advanced",
                                radioButtons("inh",
                                             "Inheritance:",
                                             choices = list("monoallelic","biallelic"),
                                             selected = "monoallelic"),
                                numericInput("prev",
                                             "Prevalence = 1 in ... (people)",
                                             min = 1,
                                             max = 1e8,
                                             value = 2500),
                                sliderInput("hetA",
                                            "Allelic heterogeneity:",
                                            min = 0,
                                            max = 1,
                                            value = 0.1),
                                sliderInput("hetG",
                                            "Genetic heterogeneity:",
                                            min = 0,
                                            max = 1,
                                            value = 1),
                                sliderInput("pen",
                                            "Penetrance:",
                                            min = 0,
                                            max = 1,
                                            value = 0.5),
                                radioButtons("CI",
                                             "Confidence:",
                                             choices = list(0.9,0.95,0.99,0.999),
                                             selected = 0.95,
                                             inline=T),
                                numericInput("popSize",
                                             "Reference population size (alleles)",
                                             min = 1,
                                             max = 1e8,
                                             value = 2*125748)
                            )
                          ),
                          br(),    
                          selectInput("race", "Population:", c("ALL"="ALL","African/African American"="AFR","Admixed American"="AMR","Ashkenazi Jewish"="ASJ","East Asian"="EAS","Finnish"="FIN","Non-Finnish European"="NFE","Other"="OTH","South Asian"="SAS"), selectize=TRUE),
                          radioButtons("format", "Download file format", choices = c("jpeg", "png", "pdf"), selected = "pdf", inline = TRUE ),
                          actionButton("goButton", "Go!"),
                          br()
                        ),
                        mainPanel(width = 9,
                                  tabsetPanel(
                                    tabPanel("GeneInfo", br(),
                                             withSpinner(DT::dataTableOutput("geneinfo",width = "100%"),type=6),
                                             br(), DT::dataTableOutput("highlight",width = "100%"),
                                             br(), DT::dataTableOutput("clinvar",width = "100%"),
                                             br(), DT::dataTableOutput("dbnsfp",width = "100%"),
                                             br(), DT::dataTableOutput("gnomad",width = "100%"),
                                             br(), DT::dataTableOutput("uniprot"),
                                             br(), DT::dataTableOutput("ccrs")
                                    ), 
                                    tabPanel("Plot", br(),
                                             downloadButton(outputId = "down", label = "Download the plot"),
                                             withSpinner(plotlyOutput("mplot", width = "100%", height = "1000px"),type=6)
                                    ),
                                    tabPanel("Debug", br(),
                                             h5("Debug Output: consurf_score() columns and structure"),
                                             verbatimTextOutput("debug_cols")
                                    )
                                  )
                        )
                    )
                ),
           tabPanel("Variant Info",
                   sidebarLayout(
                      sidebarPanel(
                        h5("This section of the plot will be updated soon!")
                      ),
                      mainPanel(
                        h5("Hello World")
                      )
                    )
           ),
           tabPanel("About",
                    includeMarkdown("About.Rmd")
           )
           )
     )
  )
)

# Add this to your UI section (if you have access to modify it)
verbatimTextOutput("debug_consurf")
