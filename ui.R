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
        padding:10px 18px; border-radius:14px; margin:12px 12px 18px 12px;
      }
      .vv-bar{display:flex;align-items:center;justify-content:space-between;gap:16px}
      .vv-logo img{height:46px; display:block; cursor:pointer; transition:opacity 0.2s;
        background:#fff; border-radius:8px; padding:4px 10px;
      }
      .vv-logo img:hover{opacity:0.85}
      .vv-nav{display:flex;gap:16px}
      .vv-nav a{color:#fff; opacity:.95; font-weight:700; padding:8px 12px; border-radius:10px; text-decoration:none}
      .vv-nav a:hover{background:#ffffff22; opacity:1}
      .vv-cta{display:inline-block; margin-top:14px; padding:10px 16px; background:var(--c5); color:#1a1a1a; border-radius:10px; font-weight:700; border:none; cursor:pointer}
      .vv-cta:hover{background:#ffb733}
      .vv-card{background:#fff; border:1px solid #d8d7d9; border-radius:16px; padding:18px; margin-bottom:18px}
      .vv-badge{display:inline-block; padding:4px 8px; border-radius:999px; font-size:12px; font-weight:700; color:#fff}
      .vv-badge.sources{background:var(--c2)} .vv-badge.outputs{background:var(--c3)}
      .vv-badge.use{background:var(--c4)} .vv-badge.quick{background:var(--c2)}
      .vv-barline{height:6px; width:72px; border-radius:6px; background:var(--c4); margin:6px 0 14px}
      #shiny-notification-panel{
        position:fixed; top:90px; right:18px; bottom:auto;
        width:280px; z-index:9999;
      }
      .shiny-notification{
        background:#fff; border:1px solid #e2e8f0;
        border-radius:10px; box-shadow:0 4px 16px rgba(0,0,0,0.12);
        padding:12px 16px; font-size:13px; color:#1e293b;
      }
      .shiny-notification::before{
        content:'Be Patient...';
        display:block; font-weight:700; font-size:13px;
        color:#FF4500; margin-bottom:4px;
      }
      .shiny-notification-close{ color:#94a3b8; }
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
  
  # header (Logo clickable to go Home)
  div(class="vv-appbar",
      div(class="vv-bar",
          div(class="vv-logo",
              tags$a(href="#", onclick="vvGo('home'); return false;",
                     tags$img(src="VarViz_Logo.png", alt="VarViz"))),
          div(class="vv-nav",
              tags$a(href="#", onclick="vvGo('protein_view'); return false;", "Plot Gene Variants"),
              tags$a(href="#", onclick="vvGo('help'); return false;", "Help")
          )
      )
  ),
  
  # ---------- HOME ----------
  conditionalPanel(condition = "input.route == 'home'",
                   div(style="padding:0 12px; max-width:1100px; margin:0 auto;",
                       
                       # --- Hero card ---
                       div(class="vv-card", style="text-align:center; padding:36px 24px 30px; border:none; background:linear-gradient(160deg,#d1dce8 0%,#bccadc 100%); position:relative; overflow:hidden;",
                           div(style="position:absolute; top:-30px; right:-30px; width:160px; height:160px; border-radius:50%; background:rgba(255,163,0,.18);"),
                           div(style="position:absolute; bottom:-20px; left:-20px; width:120px; height:120px; border-radius:50%; background:rgba(188,80,144,.16);"),
                           h1(style="font-size:36px; font-weight:800; margin:0 0 8px; letter-spacing:.3px;",
                              HTML("<span style='color:#d94030;'>Var</span><span style='color:#3d5a6e;'>Viz</span>")),
                           p(style="font-size:18px; color:#041836; margin:0 0 4px; max-width:640px; margin-left:auto; margin-right:auto;",
                             "Protein-centric gene variant visualization with integrated ACMG/AMP pathogenicity classification."),
                           p(style="font-size:14px; color:#182840; margin:0 0 20px;",
                             "Align variants to domains, structure, population frequency, 20+ computational predictors, and clinical evidence. Hybrid Richards 2015 rule-based + Tavtigian 2020 Bayesian ACMG scoring with Pejaver 2022 PP3/BP4 calibration. All data fetched live from APIs; results cached in-session."),
                           tags$button(class="vv-cta", onclick="vvGo('protein_view');",
                                       style="font-size:16px; padding:12px 28px; border-radius:12px;",
                                       HTML("&#9654;&ensp;Start Analyzing"))
                       ),
                       
                       # --- How it works (numbered steps) ---
                       div(style="margin:28px 0 8px;",
                           h2(style="color:#003f5c; font-size:22px; margin:0 0 4px;", "How It Works"),
                           div(class="vv-barline")
                       ),
                       fluidRow(
                         column(3,
                                div(class="vv-card", style="text-align:center; padding:22px 14px;",
                                    div(style="width:44px; height:44px; border-radius:12px; background:#003f5c; color:#fff; display:inline-flex; align-items:center; justify-content:center; font-weight:800; font-size:20px; margin-bottom:10px;", "1"),
                                    h4(style="margin:0 0 4px; color:#003f5c; font-size:15px;", "Select Gene"),
                                    p(style="font-size:13px; color:#64748b; margin:0;", "Search by HGNC symbol. VarViz maps it to UniProt, coordinates, and protein length.")
                                )
                         ),
                         column(3,
                                div(class="vv-card", style="text-align:center; padding:22px 14px;",
                                    div(style="width:44px; height:44px; border-radius:12px; background:#58508d; color:#fff; display:inline-flex; align-items:center; justify-content:center; font-weight:800; font-size:20px; margin-bottom:10px;", "2"),
                                    h4(style="margin:0 0 4px; color:#58508d; font-size:15px;", "Enter Variants"),
                                    p(style="font-size:13px; color:#64748b; margin:0;", "Paste in p-notation (e.g. p.R175H) or upload a file with custom positions.")
                                )
                         ),
                         column(3,
                                div(class="vv-card", style="text-align:center; padding:22px 14px;",
                                    div(style="width:44px; height:44px; border-radius:12px; background:#bc5090; color:#fff; display:inline-flex; align-items:center; justify-content:center; font-weight:800; font-size:20px; margin-bottom:10px;", "3"),
                                    h4(style="margin:0 0 4px; color:#bc5090; font-size:15px;", "Visualize"),
                                    p(style="font-size:13px; color:#64748b; margin:0;", "Click Go. VarViz fetches from 9+ live APIs (UniProt, AlphaFold, gnomAD, ClinVar, CCRS, dbNSFP, UCSC, Ensembl) with progress updates. Cached genes load without re-fetching.")
                                )
                         ),
                         column(3,
                                div(class="vv-card", style="text-align:center; padding:22px 14px;",
                                    div(style="width:44px; height:44px; border-radius:12px; background:#ff6361; color:#fff; display:inline-flex; align-items:center; justify-content:center; font-weight:800; font-size:20px; margin-bottom:10px;", "4"),
                                    h4(style="margin:0 0 4px; color:#ff6361; font-size:15px;", "Interpret & Export"),
                                    p(style="font-size:13px; color:#64748b; margin:0;", "Explore GeneInfo tab, hover interactive plots, and export publication-ready PDF, PNG, or JPEG.")
                                )
                         )
                       ),
                       
                       # --- Feature cards (2 cols) ---
                       div(style="margin:24px 0 8px;",
                           h2(style="color:#003f5c; font-size:22px; margin:0 0 4px;", "What You Get"),
                           div(class="vv-barline")
                       ),
                       fluidRow(
                         column(6,
                                div(class="vv-card",
                                    div(style="display:flex; align-items:center; gap:10px; margin-bottom:10px;",
                                        div(style="width:36px; height:36px; border-radius:10px; background:linear-gradient(135deg,#003f5c,#58508d); display:flex; align-items:center; justify-content:center;",
                                            HTML("<span style='color:#fff; font-size:18px;'>&#x1F4CA;</span>")),
                                        h3(style="margin:0; font-size:17px; color:#003f5c;", "9 Interactive Tracks")
                                    ),
                                    tags$ul(style="margin:0; padding-left:20px; color:#041836; font-size:14px;",
                                      tags$li("AlphaFold pLDDT structural confidence"),
                                      tags$li("gnomAD allele frequency rainfall plot"),
                                      tags$li("Mutation density (gnomAD + ClinVar + user variants)"),
                                      tags$li("AlphaFold Mean Pathogenicity scores"),
                                      tags$li("ClinVar / PTM / CCRS combined track"),
                                      tags$li(HTML("<strong>Multi Conservation (UCSC)</strong>: PhyloP 100V, PhyloP 470M, PhastCons via UCSC REST API")),
                                      tags$li("Conservation scores (ConSurf upload)"),
                                      tags$li("UniProt domains, regions & features"),
                                      tags$li("Transmembrane helices & topological domains")
                                    )
                                )
                         ),
                         column(6,
                                div(class="vv-card",
                                    div(style="display:flex; align-items:center; gap:10px; margin-bottom:10px;",
                                        div(style="width:36px; height:36px; border-radius:10px; background:linear-gradient(135deg,#bc5090,#ff6361); display:flex; align-items:center; justify-content:center;",
                                            HTML("<span style='color:#fff; font-size:18px;'>&#x1F50D;</span>")),
                                        h3(style="margin:0; font-size:17px; color:#bc5090;", "Rich Gene Info")
                                    ),
                                    tags$ul(style="margin:0; padding-left:20px; color:#041836; font-size:14px;",
                                      tags$li("Variant Annotation Summary table with 40+ columns across all tracks"),
                                      tags$li(HTML("dbNSFP pathogenicity scores via MyVariant.info - SIFT, PolyPhen2, REVEL, CADD, FATHMM, PROVEAN, MutationTaster & more")),
                                      tags$li(HTML("ACMG evidence tags auto-computed: hybrid Richards 2015 rule-based + Tavtigian 2020 Bayesian scoring + Pejaver 2022 PP3/BP4 calibration")),
                                      tags$li(HTML("Conservation: ScoreCons (CCRStoAAC), GERP++, PhyloP, PhastCons")),
                                      tags$li("Population frequencies by ancestry (gnomAD, 1000G, ExAC)"),
                                      tags$li("AlphaMissense exact per-variant scores from AlphaFold substitution matrix"),
                                      tags$li("Two-tier ClinVar matching: exact variant vs. same-position context"),
                                      tags$li("Gene summary card, disease associations, word cloud"),
                                      tags$li("15 external database links (OMIM, ClinGen, DECIPHER, GTEx, ConSurf-DB, PanelApp...)"),
                                      tags$li("In-session caching: genes re-load without re-fetching")
                                    )
                                )
                         )
                       ),
                       
                       # --- Data sources + Use cases side by side ---
                       fluidRow(
                         column(6,
                                div(class="vv-card",
                                    div(style="display:flex; align-items:center; gap:10px; margin-bottom:10px;",
                                        div(style="width:36px; height:36px; border-radius:10px; background:linear-gradient(135deg,#58508d,#bc5090); display:flex; align-items:center; justify-content:center;",
                                            HTML("<span style='color:#fff; font-size:18px;'>&#x1F5C4;</span>")),
                                        h3(style="margin:0; font-size:17px; color:#58508d;", "Data Sources")
                                    ),
                                    tags$ul(style="margin:0; padding-left:20px; color:5569; font-size:14px;",
                                      tags$li(HTML("<strong>UniProt REST API</strong> - domains, features, PTMs, function, disease")),
                                      tags$li(HTML("<strong>AlphaFold API</strong> - per-residue pLDDT & mean pathogenicity; AlphaMissense substitution matrix")),
                                      tags$li(HTML("<strong>gnomAD GraphQL API</strong> - population frequencies (v4 + v2 fallback)")),
                                      tags$li(HTML("<strong>NCBI E-utilities</strong> - ClinVar pathogenic variants & gene summaries")),
                                      tags$li(HTML("<strong>UCSC REST API</strong> - PhyloP 100V, PhyloP 470M, PhastCons (hg38) for Multi Conservation track")),
                                      tags$li(HTML("<strong>Ensembl REST API</strong> - canonical transcript exon structure for UCSC score mapping")),
                                      tags$li(HTML("<strong>GENCODE v49</strong> - gene coordinates (GRCh38)")),
                                      tags$li(HTML("<strong>CCRS</strong> - constrained coding regions (local indexed file)"))
                                    )
                                )
                         ),
                         column(6,
                                div(class="vv-card",
                                    div(style="display:flex; align-items:center; gap:10px; margin-bottom:10px;",
                                        div(style="width:36px; height:36px; border-radius:10px; background:linear-gradient(135deg,#ff6361,#ffa600); display:flex; align-items:center; justify-content:center;",
                                            HTML("<span style='color:#fff; font-size:18px;'>&#x1F3AF;</span>")),
                                        h3(style="margin:0; font-size:17px; color:#ff6361;", "Use Cases")
                                    ),
                                    tags$ul(style="margin:0; padding-left:20px; color:#041836; font-size:14px;",
                                      tags$li("Triage candidate variants for clinical review"),
                                      tags$li("Identify pathogenic hotspots and benign deserts"),
                                      tags$li("Apply allele frequency filters based on disease architecture"),
                                      tags$li("Create publication-ready protein landscape figures"),
                                      tags$li("Investigate gene-disease associations with external links"),
                                      tags$li("Variant selection for mutational experiments, particularly Deep Mutational Scanning (DMS)")
                                    )
                                )
                         )
                       ),
                       
                       # --- Bottom CTA ---
                       div(style="text-align:center; margin:24px 0 12px;",
                           tags$button(class="vv-cta", onclick="vvGo('protein_view');",
                                       style="font-size:15px; padding:11px 24px; border-radius:12px; margin-right:12px;",
                                       HTML("&#9654;&ensp;Plot Gene Variants")),
                           tags$button(class="vv-cta", onclick="vvGo('help');",
                                       style="font-size:15px; padding:11px 24px; border-radius:12px; background:#58508d; color:#fff;",
                                       HTML("&#x1F4D6;&ensp;Read Tutorial"))
                       ),
                       
                       # --- Footer ---
                       div(style="margin-top:20px; padding:14px 0; border-top:1px solid #e2e8f0; text-align:center; color:#182840; font-size:13px;",
                           HTML("VarViz - Protein-centric gene variant visualization - 2026")
                       )
                   )
  ),
  
  # ---------- PLOT GENE VARIANTS (was protein_view) ----------
  conditionalPanel(condition = "input.route == 'protein_view'",
                   fluidRow(
                     column(3,
                            tags$div(style = "background:#ffffff; border:1px solid #e2e8f0; border-radius:14px; padding:16px 14px 10px; box-shadow:0 1px 4px rgba(0,0,0,0.04);",
                            uiOutput("selects1"),
                            selectizeInput("gene_name", "Choose a gene", choices = NULL),
                            fileInput("variant_file",
                                      label = NULL,
                                      buttonLabel = "Upload variants (.txt)",
                                      placeholder = "or paste below",
                                      accept = c("text/plain", ".txt", ".csv")
                            ),
                            textAreaInput(
                              "variants",
                              "Variants (comma-separated, p-notation):",
                              "p.I554N",
                              height = "100px"
                            ),
                            tags$div(style = "text-align:center; margin:10px 0 6px;",
                              actionButton("goButton", "Go!",
                                           style = "background:#003f5c; color:#fff; font-weight:700; width:75%; max-width:200px; border-radius:8px; padding:8px 0; font-size:15px;")
                            ),
                            tags$div(
                              style = "text-align:center; margin:0 0 14px;",
                              tags$span(
                                HTML("&#9432; Gene &amp; variant(s) are required"),
                                style = paste0(
                                  "font-size:11px; color:#94a3b8; font-style:italic;",
                                  "display:inline-flex; align-items:center; gap:4px;"
                                )
                              )
                            ),
                            checkboxGroupInput(
                              "plotselection", "Output Tracks:",
                              choices = c("gnomAD Frequency"          = "freq",
                                          "Mutation Density"          = "density",
                                          "AF Mean Pathogenicity"     = "clinvar",
                                          "ClinVar / PTM / CCRS"      = "clinvar_ptm",
                                          "Multi Conservation (UCSC)" = "multiconservation"),
                              selected = c("freq","density","clinvar","clinvar_ptm","multiconservation")
                            ),
                            # --- Frequency Cutoff Method ---
                            selectInput("cutoff_method", 
                                        tags$span(style = "font-weight:700; color:#003f5c;", 
                                                  "Frequency Cutoff Method:"),
                                        choices = list(
                                          "Calculate AF (by prevalence)" = "calc_af",
                                          "AC Filter (by max population AF)" = "ac_filter"
                                        ),
                                        selected = "calc_af"),
                            # Result summary (always visible)
                            htmlOutput("cutoff_summary"),
                            # Toggle link
                            tags$a(id = "toggleCutoff", 
                                   style = "cursor:pointer; font-size:12px; color:#58508d; font-weight:600;",
                                   icon("sliders-h"), " Hide parameters"),
                            br(),
                            # Parameter detail panels — shown by default, toggle to hide
                            div(id = "cutoff_details",
                                  # --- Calculate AF Panel ---
                                  conditionalPanel(
                                    condition = "input.cutoff_method == 'calc_af'",
                                    tags$div(
                                      style = "background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:14px; margin-top:8px; margin-bottom:10px;",
                                      tags$p(style = "font-size:11px; color:#64748b; margin:0 0 8px;",
                                             tags$a(href = "https://www.sciencedirect.com/science/article/pii/S1098360021013678",
                                                    target = "_blank", "Ware et al. (2018) Genetics in Medicine")),
                                      radioButtons("inh","Inheritance:", inline = TRUE,
                                                   choices = list("monoallelic","biallelic"),
                                                   selected = "monoallelic"),
                                      numericInput("prev","Prevalence = 1 in ... (people)",
                                                   min = 1, max = 1e8, value = 2000),
                                      sliderInput("hetA","Allelic heterogeneity:",
                                                  min = 0, max = 1, value = 0.2, step = 0.01),
                                      sliderInput("hetG","Genetic heterogeneity:",
                                                  min = 0, max = 1, value = 1, step = 0.01),
                                      sliderInput("pen","Penetrance:",
                                                  min = 0, max = 1, value = 0.5, step = 0.01),
                                      selectInput("CI_af", "Confidence:",
                                                  choices = list("0.9" = 0.9, "0.95" = 0.95,
                                                                 "0.99" = 0.99, "0.999" = 0.999),
                                                  selected = 0.95),
                                      numericInput("popSize_af","Reference population size (alleles):",
                                                   min = 1, max = 1e8, value = 2*125748)
                                    )
                                  ),
                                  # --- AC Filter Panel ---
                                  conditionalPanel(
                                    condition = "input.cutoff_method == 'ac_filter'",
                                    tags$div(
                                      style = "background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:14px; margin-top:8px; margin-bottom:10px;",
                                      tags$p(style = "font-size:11px; color:#64748b; margin:0 0 8px;",
                                             tags$a(href = "https://cardiodb.org/allelefrequencyapp/",
                                                    target = "_blank", "CardioDB Allele Frequency App")),
                                      numericInput("maxPopAF", "Maximum population AF:",
                                                   min = 0, max = 1, value = 0.001, step = 0.0001),
                                      numericInput("popSize","Reference population size (alleles):",
                                                   min = 1, max = 1e8, value = 2*125748),
                                      selectInput("CI", "Confidence:",
                                                  choices = list("0.9" = 0.9, "0.95" = 0.95,
                                                                 "0.99" = 0.99, "0.999" = 0.999),
                                                  selected = 0.95)
                                    )
                                  )
                            ),
                            br(),
                            fileInput("consurf_file",
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".txt"),
                                      label = h5("ConSurf grades file (tab-separated):")
                            ),
                            # Custom pathogenic positions file removed — redundant with variant input
                            radioButtons("label", "Variant Label:", inline = TRUE,
                                         choices = list("Yes","No"), selected = "Yes"),
                            radioButtons("clinvar_filter", "ClinVar Variants Filter:", inline = TRUE,
                                         choices = list("Any Star","1 Star or more","2 Stars or more"),
                                         selected = "1 Star or more"),
                            # De Novo Status: moved to per-variant card dropdowns
                            # (applied individually per variant, not globally)
                            radioButtons("format", "Download format:",
                                         choices = c("jpeg","png","pdf"), selected = "pdf", inline = TRUE),
                            br()
                            )
                     ),
                     column(9,
                            tabsetPanel(id = "main_tabs",
                              tabPanel("GeneInfo", br(),
                                       # Top row: Gene summary (left) + Word Cloud (right)
                                       fluidRow(
                                         column(6, htmlOutput("geneinfo_summary")),
                                         column(6,
                                                div(style="background:#fff; border:1px solid #d8d7d9; border-radius:16px; padding:18px; margin-bottom:18px; height:100%;",
                                                    div(style="display:flex; align-items:center; margin-bottom:12px;",
                                                        div(style="height:6px; width:72px; border-radius:6px; background:#58508d; margin-right:14px;"),
                                                        h4(style="color:#003f5c; margin:0; font-size:16px;", "Gene Word Cloud")
                                                    ),
                                                    p(style="font-size:12px; color:#182840; margin:0 0 8px;",
                                                      "Keywords from UniProt & NCBI annotations"),
                                                    plotOutput("wordcloud_plot", height = "280px")
                                                )
                                         )
                                       ),
                                       # Remaining cards: Function, Disease, Links
                                       htmlOutput("geneinfo_details"),
                                       # Variant Intersection Table (appears after data loads)
                                       uiOutput("variant_table_section"),
                                       # Download button for variant table TSV
                                       conditionalPanel(
                                         condition = "output.variant_table_ready",
                                         div(style = "margin:-10px 0 18px; text-align:right;",
                                             downloadButton("download_gnomad_raw", 
                                                            label = "Download gnomAD data",
                                                            style = "margin-right:8px; background:#1d6a96; color:#fff; border:none; border-radius:6px; padding:5px 12px; font-size:12px;"),
                                             downloadButton("download_variant_table", 
                                                            label = "Download TSV",
                                                            style = "background:#003f5c; color:#fff; border:none; border-radius:8px; padding:8px 18px; font-size:13px; cursor:pointer;")
                                         )
                                       ),
                                       uiOutput("external_links_section")
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
                       tags$iframe(src = "help.html", width = "100%", height = "900px",
                                  style = "border:none; border-radius:12px;")
                   )
  )
))
