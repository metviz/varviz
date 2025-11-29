# ui.R
suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinycssloaders))
suppressMessages(library(DT))
suppressMessages(library(plotly))
suppressMessages(library(markdown))

shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "shortcut icon", href = "favicon.ico"),
    tags$style(HTML("
      :root{ --c1:#003f5c; --c2:#58508d; --c3:#bc5090; --c4:#ff6361; --c5:#ffa600; --bg:#f9fafb; --ink:#0f172a; }
      html,body{background:var(--bg); color:var(--ink)}
      .vv-appbar{
        background:linear-gradient(135deg,var(--c1),var(--c2));
        color:#fff; border-bottom:6px solid var(--c5);
        padding:14px 18px; border-radius:14px; margin:12px 12px 18px 12px;
      }
      .vv-bar{display:flex;align-items:center;justify-content:space-between;gap:16px}
      .vv-title{font-weight:800;font-size:28px;letter-spacing:.3px;margin:0}
      .vv-title a{color:#fff;text-decoration:none}
      .vv-title a:hover{opacity:.9}
      .vv-nav{display:flex;gap:16px}
      .vv-nav a{color:#fff; opacity:.95; font-weight:700; padding:8px 12px; border-radius:10px; text-decoration:none}
      .vv-nav a:hover{background:#ffffff22; opacity:1}
      .vv-cta{display:inline-block; margin-top:14px; padding:10px 16px; background:var(--c5); color:#1a1a1a; border-radius:10px; font-weight:700; border:none; cursor:pointer}
      .vv-cta:hover{background:#ffb733}
      .vv-card{background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:18px; margin-bottom:18px}
      .vv-badge{display:inline-block; padding:4px 8px; border-radius:999px; font-size:12px; font-weight:700; color:#fff}
      .vv-badge.sources{background:var(--c2)} .vv-badge.outputs{background:var(--c3)}
      .vv-badge.use{background:var(--c4)} .vv-badge.quick{background:var(--c2)}
      .vv-barline{height:6px; width:72px; border-radius:6px; background:var(--c4); margin:6px 0 14px}
    ")),
    tags$script(HTML("
      window.vvGo = function(route){
        Shiny.setInputValue('route', route, {priority:'event'});
        window.scrollTo({top:0, behavior:'smooth'});
      };
      $(function(){
        if(!Shiny.shinyapp || !Shiny.shinyapp.$inputValues || !Shiny.shinyapp.$inputValues.route){
          Shiny.setInputValue('route', 'home', {priority:'event'});
        }
      });
    "))
  ),
  
  useShinyjs(),
  
  # hidden route value so conditionalPanels work at load
  tags$div(style="display:none", textInput("route", label = NULL, value = "home")),
  
  # header (VarViz clickable to go Home)
  div(class="vv-appbar",
      div(class="vv-bar",
          tags$h1(class="vv-title",
                  tags$a(href="#", onclick="vvGo('home'); return false;", "VarViz")),
          div(class="vv-nav",
              tags$a(href="#", onclick="vvGo('protein_view'); return false;", "Plot Gene Variants"),
              tags$a(href="#", onclick="vvGo('help'); return false;", "Help")
          )
      )
  ),
  
  # ---------- HOME ----------
  conditionalPanel(condition = "input.route == 'home'",
                   div(style="padding:0 12px",
                       div(class="vv-card",
                           h2(style="color:#003f5c;margin:0 0 6px","VarViz"),
                           p("Protein-centric variant visualization for fast, informed calls."),
                           # actionButton("launchApp", "Launch Plot Gene Variants", class = "vv-cta",
                           #              onclick = "vvGo('protein_view');")
                       ),
                       div(class="vv-barline"),
                       h2(style="color:#003f5c;margin-top:0","Overview"),
                       p("VarViz aligns variants to protein features so you can spot hotspots, deserts, and patterns in seconds. Built for quick review, clean exports, and team sharing."),
                       
                       fluidRow(
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge", style="background:#ffa600;color:#1a1a1a","What You Can Do"),
                                    tags$ul(
                                      tags$li("Map variants on a shared residue axis with domains and motifs."),
                                      tags$li("Review structure confidence with AlphaFold pLDDT."),
                                      tags$li("Check conservation and population frequency."),
                                      tags$li("Scan ClinVar and dbNSFP pathogenicity trends."),
                                      tags$li("Export publication-ready figures.")
                                    )
                                )
                         ),
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge sources","Data Sources"),
                                    tags$ul(
                                      tags$li("UniProt and Pfam for protein features."),
                                      tags$li("AlphaFold pLDDT for model confidence."),
                                      tags$li("gnomAD for frequency and density."),
                                      tags$li("ClinVar summaries for clinical assertions."),
                                      tags$li("dbNSFP, including REVEL, for pathogenicity scores.")
                                    )
                                )
                         )
                       ),
                       
                       fluidRow(
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge", style="background:#003f5c","How It Works"),
                                    tags$ol(
                                      tags$li("Choose a gene or paste variants."),
                                      tags$li("Fetch features and scores from APIs or local files."),
                                      tags$li("Render aligned tracks on a shared residue axis."),
                                      tags$li("Interact, filter, annotate, and export.")
                                    )
                                )
                         ),
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge outputs","Outputs"),
                                    tags$ul(
                                      tags$li("Interactive stacked tracks."),
                                      tags$li("PNG and PDF figure exports."),
                                      tags$li("Optional pLDDT JSON downloads."),
                                      tags$li("Consistent captions and styling.")
                                    ),
                                    div(style="margin-top:12px",
                                        # actionButton("launchApp2", "Open Plot Gene Variants", class = "vv-cta",
                                        #              onclick = "vvGo('protein_view');")
                                    )
                                )
                         )
                       ),
                       
                       fluidRow(
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge use","Use Cases"),
                                    tags$ul(
                                      tags$li("Triage variants for tumor boards."),
                                      tags$li("Highlight benign deserts and pathogenic hotspots."),
                                      tags$li("Create figures for manuscripts and posters."),
                                      tags$li("Teach protein context in genetics courses.")
                                    )
                                )
                         ),
                         column(6,
                                div(class="vv-card",
                                    span(class="vv-badge quick","Quick Start"),
                                    tags$ol(
                                      tags$li("Open the app."),
                                      tags$li("Select a gene."),
                                      tags$li("Paste variants or upload a file."),
                                      tags$li("Pick tracks and export.")
                                    )
                                )
                         )
                       )
                   )
  ),
  
  # ---------- PLOT GENE VARIANTS (was protein_view) ----------
  conditionalPanel(condition = "input.route == 'protein_view'",
                   fluidRow(
                     column(3,
                            uiOutput("selects1"),
                            checkboxGroupInput(
                              "plotselection", "Output Plot:",
                              choices = c("Clinvar" = "clinvar",
                                          "Density" = "density",
                                          "gnomAD Freq" = "freq",
                                          "gnomAD Depth" = "depth"),
                              selected = c("clinvar","density","freq","depth")
                            ),
                            selectizeInput("gene_name", "Choose a gene", choices = NULL),
                            fileInput("consurf_file",
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".txt"),
                                      label = h5("Please provide Consurf scores in tab separated format:")
                            ),
                            fileInput("user_file",
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".txt"),
                                      label = h5("Please provide protein locations:")
                            ),
                            textAreaInput(
                              "variants",
                              "Please provide comma separated variants in pMut format:",
                              "p.I554N",
                              height = "100px"
                            ),
                            radioButtons("label", "Variant Label:", inline = TRUE,
                                         choices = list("Yes","No"), selected = "Yes"),
                            radioButtons("clinvar_filter", "Clinvar Variants Filter :", inline = TRUE,
                                         choices = list("Any Star","1 Star or more","2 Stars or more"),
                                         selected = "1 Star or more"),
                            htmlOutput("maxAC_F"),
                            br(),
                            a(id = "toggleAdvanced", "Show/Hide this section", href = "#"), br(),
                            shinyjs::hidden(
                              div(id = "advanced",
                                  radioButtons("inh","Inheritance:",
                                               choices = list("monoallelic","biallelic"),
                                               selected = "monoallelic"),
                                  numericInput("prev","Prevalence = 1 in ... (people)", min = 1, max = 1e8, value = 2500),
                                  sliderInput("hetA","Allelic heterogeneity:", min = 0, max = 1, value = 0.1),
                                  sliderInput("hetG","Genetic heterogeneity:", min = 0, max = 1, value = 1),
                                  sliderInput("pen","Penetrance:", min = 0, max = 1, value = 0.5),
                                  radioButtons("CI","Confidence:", choices = list(0.9,0.95,0.99,0.999),
                                               selected = 0.95, inline=TRUE),
                                  numericInput("popSize","Reference population size (alleles)",
                                               min = 1, max = 1e8, value = 2*125748)
                              )
                            ),
                            br(),
                            selectInput(
                              "race", "Population:",
                              c("ALL"="ALL","African/African American"="AFR","Admixed American"="AMR",
                                "Ashkenazi Jewish"="ASJ","East Asian"="EAS","Finnish"="FIN",
                                "Non-Finnish European"="NFE","Other"="OTH","South Asian"="SAS"),
                              selectize = TRUE
                            ),
                            radioButtons("format", "Download file format",
                                         choices = c("jpeg","png","pdf"), selected = "pdf", inline = TRUE),
                            actionButton("goButton", "Go!"),
                            br()
                     ),
                     column(9,
                            tabsetPanel(
                              tabPanel("GeneInfo", br(),
                                       withSpinner(DT::dataTableOutput("geneinfo", width = "100%"), type = 6),
                                       br(), DT::dataTableOutput("highlight", width = "100%"),
                                       br(), DT::dataTableOutput("clinvar", width = "100%"),
                                       br(), DT::dataTableOutput("dbnsfp", width = "100%"),
                                       br(), DT::dataTableOutput("gnomad", width = "100%"),
                                       br(), DT::dataTableOutput("uniprot"),
                                       br(), DT::dataTableOutput("ccrs")
                              ),
                              tabPanel("Plot", br(),
                                       downloadButton(outputId = "down", label = "Download the plot"),
                                       withSpinner(plotlyOutput("mplot", width = "100%", height = "1000px"), type = 6)
                              )
                            )
                     )
                   )
  ),
  
  # ---------- HELP ----------
  conditionalPanel(condition = "input.route == 'help'",
                   div(style="padding:0 12px",
                       includeHTML("help.html")
                   )
  )
))
