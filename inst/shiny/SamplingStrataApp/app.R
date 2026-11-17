# Ensure a CRAN mirror is always set — required on shinyapps.io when any
# package function calls require() / install.packages() internally.
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ─── Package availability check ───────────────────────────────────────────────
# Missing packages are installed from CRAN before being loaded. The explicit
# library() calls below are kept (instead of a loop) so that rsconnect can still
# detect the dependencies when deploying to shinyapps.io.
required_pkgs <- c("shiny", "SamplingStrata", "DT", "ggplot2",
                   "corrplot", "shinydashboard", "shinyjs")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  message("Missing package(s): ", paste(missing_pkgs, collapse = ", "),
          " - attempting installation from CRAN...")
  tryCatch(
    install.packages(missing_pkgs, quiet = TRUE),
    error   = function(e) message("Installation error: ", conditionMessage(e)),
    warning = function(w) message("Installation warning: ", conditionMessage(w))
  )
  still_missing <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop("The following required package(s) could not be installed: ",
         paste(still_missing, collapse = ", "),
         ". Please install them manually with install.packages().",
         call. = FALSE)
  }
}

library(shiny)
library(SamplingStrata)
library(DT)
library(ggplot2)
library(corrplot)
library(shinydashboard)
library(shinyjs)

# ─── UI ────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "yellow",
  dashboardHeader(title = "SamplingStrata Workflow"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Data",            tabName = "tab_data",        icon = icon("database")),
      menuItem("2. Variables",       tabName = "tab_vars",        icon = icon("sliders-h")),
      menuItem("3. Frame & Strata",  tabName = "tab_frame",       icon = icon("layer-group")),
      menuItem("4. Constraints",     tabName = "tab_errors",      icon = icon("bullseye")),
      menuItem("5. Take-all Strata", tabName = "tab_cens",        icon = icon("filter")),
      menuItem("6. Optimization",    tabName = "tab_optim",       icon = icon("cogs")),
      menuItem("7. Results",         tabName = "tab_results",     icon = icon("chart-bar")),
      menuItem("8. Sample",          tabName = "tab_sample",      icon = icon("table")),
      menuItem("9. R Script",        tabName = "tab_script",      icon = icon("file-code"))
    ),
    div(style = "padding:12px 15px; border-top: 1px solid #444; margin-top:10px;",
      actionButton("btn_help", " Help Guide", icon = icon("question-circle"),
                   class = "btn-info btn-block",
                   style = "font-weight:bold;")
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .nav-tabs-custom .nav-tabs li.active { border-top-color: #605ca8; }
      .box-body .dataTables_wrapper { overflow-x: auto; }
      pre { max-height: 300px; overflow-y: auto; }
    "))),
    tabItems(

      # ── TAB 1: DATA ──────────────────────────────────────────────────────────
      tabItem(tabName = "tab_data",
        fluidRow(
          box(title = "Data Source", width = 4, status = "primary", solidHeader = TRUE,
            radioButtons("data_source", "Choose data:",
              choices = c("Use swissmunicipalities (demo)" = "demo",
                          "Upload CSV file"               = "upload"),
              selected = "demo"),
            conditionalPanel("input.data_source == 'upload'",
              fileInput("upload_file", "Upload CSV", accept = ".csv"),
              checkboxInput("header", "Header row", TRUE),
              radioButtons("sep", "Separator",
                choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), inline = TRUE)
            ),
            actionButton("btn_load", "Load Data", class = "btn-primary btn-block",
                         icon = icon("upload"))
          ),
          box(title = "Dataset Preview", width = 8, status = "info",
            DTOutput("tbl_data")
          )
        ),
        fluidRow(
          box(title = "Dataset Summary", width = 12,
            verbatimTextOutput("txt_summary")
          )
        )
      ),

      # ── TAB 2: VARIABLES ─────────────────────────────────────────────────────
      tabItem(tabName = "tab_vars",
        fluidRow(
          box(title = "Variable Assignment", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("ui_id_var"),
            uiOutput("ui_domain_var"),
            uiOutput("ui_y_vars"),
            uiOutput("ui_x_vars"),
            uiOutput("ui_x_bins_table"),
            actionButton("btn_prepare", "Prepare Variables", class = "btn-success btn-block",
                         icon = icon("check"))
          ),
          box(title = "Correlation Matrix (X vs Y)", width = 8, status = "info",
            plotOutput("plot_corr", height = "420px")
          )
        ),
        fluidRow(
          box(title = "X variables after categorisation (head)", width = 12,
            DTOutput("tbl_xbins")
          )
        ),
        fluidRow(
          box(title = "Categorised X variables — frequency table", width = 12,
            uiOutput("ui_xfreq_note"),
            DTOutput("tbl_xfreq")
          )
        )
      ),

      # ── TAB 3: FRAME & STRATA ────────────────────────────────────────────────
      tabItem(tabName = "tab_frame",
        fluidRow(
          box(title = "Build Frame & Strata", width = 3, status = "primary", solidHeader = TRUE,
            p("Uses the variable assignments from the previous step."),
            actionButton("btn_build_frame", "Build Frame", class = "btn-primary btn-block",
                         icon = icon("play")),
            br(),
            actionButton("btn_build_strata", "Build Strata", class = "btn-success btn-block",
                         icon = icon("play"))
          ),
          box(title = "Frame (head)", width = 9, status = "info",
            DTOutput("tbl_frame")
          )
        ),
        fluidRow(
          box(title = "Strata (head)", width = 12, status = "info",
            DTOutput("tbl_strata")
          )
        )
      ),

      # ── TAB 4: CONSTRAINTS ───────────────────────────────────────────────────
      tabItem(tabName = "tab_errors",
        fluidRow(
          # ── LEFT: global defaults ──────────────────────────────────────────
          box(title = "Default CV (all domains)", width = 3,
              status = "primary", solidHeader = TRUE,
            p("Set a default CV for each target variable, then apply it to all domains at once."),
            uiOutput("ui_cv_global"),
            actionButton("btn_apply_global_cv", "Apply to all domains",
                         class = "btn-primary btn-block", icon = icon("arrow-right")),
            hr(),
            actionButton("btn_set_errors", "Set Constraints", class = "btn-success btn-block",
                         icon = icon("check")),
            br(),
            actionButton("btn_bethel", "Bethel Baseline", class = "btn-info btn-block",
                         icon = icon("calculator")),
            br(),
            actionButton("btn_check", "Check Input", class = "btn-warning btn-block",
                         icon = icon("stethoscope")),
            hr(),
            h5(strong("Import / Export")),
            downloadButton("btn_download_cv", "Download CV as CSV",
                           class = "btn-default btn-block", icon = icon("download")),
            br(),
            fileInput("upload_cv", "Upload CV CSV", accept = ".csv",
                      buttonLabel = icon("upload"), placeholder = "No file selected"),
            uiOutput("ui_cv_upload_status")
          ),
          # ── RIGHT: per-domain detail + results ────────────────────────────
          box(title = "Per-domain CV constraints", width = 9,
              status = "info", solidHeader = TRUE,
            p("Fine-tune the CV for each variable and domain. Click \"Apply to all domains\" on the left to reset all values to the global defaults."),
            uiOutput("ui_errors_form"),
            hr(),
            h4("Errors Dataframe"),
            DTOutput("tbl_errors"),
            br(),
            h4("Bethel Baseline Total Sample Size"),
            verbatimTextOutput("txt_bethel")
          )
        ),
        fluidRow(
          box(title = "Check Input Result", width = 12,
            verbatimTextOutput("txt_check")
          )
        )
      ),

      # ── TAB 5: TAKE-ALL STRATA ───────────────────────────────────────────────
      tabItem(tabName = "tab_cens",
        fluidRow(
          box(title = "Take-all Strata Settings (optional)", width = 4,
              status = "primary", solidHeader = TRUE,
            p("Optional step. Upload a CSV file containing the units that must be",
              "included with certainty. The file must be a subset of the dataset",
              "loaded in Step 1."),
            checkboxInput("use_cens", "Enable take-all strata", value = FALSE),
            conditionalPanel("input.use_cens == true",
              hr(),
              fileInput("cens_file", "Upload take-all subset CSV", accept = ".csv"),
              checkboxInput("cens_header", "Header row", TRUE),
              radioButtons("cens_sep", "Separator",
                choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), inline = TRUE),
              actionButton("btn_apply_cens", "Validate and build take-all frame",
                           class = "btn-warning btn-block", icon = icon("upload"))
            )
          ),
          box(title = "Take-all units — summary & preview", width = 8,
              status = "info",
            verbatimTextOutput("txt_cens_preview"),
            br(),
            DTOutput("tbl_cens_head")
          )
        )
      ),

      # ── TAB 6: OPTIMIZATION ─────────────────────────────────────────────────
      tabItem(tabName = "tab_optim",
        fluidRow(
          box(title = "Optimisation Settings", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("ui_optim_method"),
            numericInput("nstrata", "Max strata per domain (maxclusters for KMeans)", 10, 2, 50),
            checkboxInput("kmeans_nstrata", "Use nStrata from KmeansSolution", FALSE),
            uiOutput("ui_kmeans_nstrata_info"),
            numericInput("iter",    "Iterations (iter)",    50,  5, 500),
            numericInput("pops",    "Population size (pops)", 10, 5, 50),
            conditionalPanel("input.optim_method == 'atomic'",
              checkboxInput("use_kmeans", "Use KmeansSolution as seed", TRUE)
            ),
            conditionalPanel("input.optim_method == 'continuous'",
              checkboxInput("use_suggestions",
                            "Use suggestions (prepareSuggestion)", FALSE)
            ),
            actionButton("btn_optim", "Run Optimisation", class = "btn-danger btn-block",
                         icon = icon("rocket")),
            br(),
            p(em("This may take a while. Output appears below when done."))
          ),
          box(title = "Optimisation Log", width = 8, status = "warning",
            verbatimTextOutput("txt_optim_log"),
            hr(),
            h4("Convergence Plot"),
            fluidRow(
              column(2, actionButton("btn_plot_prev", "", icon = icon("chevron-left"),
                                    class = "btn-sm btn-default")),
              column(8, uiOutput("txt_plot_page", inline = TRUE),
                     style = "text-align:center; padding-top:6px;"),
              column(2, actionButton("btn_plot_next", "", icon = icon("chevron-right"),
                                    class = "btn-sm btn-default"))
            ),
            div(style = "width:100%; margin-bottom:10px;",
              imageOutput("plot_optim", width = "100%", height = "auto")
            ),
            hr(),
            h4("Expected CV after optimisation"),
            verbatimTextOutput("txt_expected_cv"),
            hr(),
            h4("Strata Summary"),
            DTOutput("tbl_strata_summary")
          )
        )
      ),

      # ── TAB 7: RESULTS ───────────────────────────────────────────────────────
      tabItem(tabName = "tab_results",
        fluidRow(
          box(title = "Evaluation Settings", width = 3, status = "primary", solidHeader = TRUE,
            numericInput("nsampl", "Simulations (nsampl)", 50, 10, 500),
            actionButton("btn_eval", "Evaluate Solution", class = "btn-primary btn-block",
                         icon = icon("flask")),
            br(),
            numericInput("adj_size", "Adjust total sample size (0 = no adjustment)",
                         0, 0, 1e6),
            actionButton("btn_adjust", "Adjust Size", class = "btn-info btn-block",
                         icon = icon("balance-scale"))
          ),
          box(title = "Strata Visualisation", width = 9, status = "info",
            uiOutput("ui_plot_vars"),
            plotOutput("plot_strata2d", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Evaluation: Coefficients of Variation", width = 6,
            plotOutput("plot_eval_cv", height = "300px")
          ),
          box(title = "Evaluation: Relative Bias", width = 6,
            plotOutput("plot_eval_bias", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Adjusted Strata (if applied)", width = 12,
            DTOutput("tbl_adj_strata")
          )
        )
      ),

      # ── TAB 8: SAMPLE ────────────────────────────────────────────────────────
      tabItem(tabName = "tab_sample",
        fluidRow(
          box(title = "Sample Selection", width = 3, status = "primary", solidHeader = TRUE,
            selectInput("sample_method", "Method",
              choices = c("Random (selectSample)" = "random",
                          "Systematic (selectSampleSystematic)" = "systematic"),
              selected = "random"),
            uiOutput("ui_sort_var"),
            actionButton("btn_select", "Select Sample", class = "btn-success btn-block",
                         icon = icon("magic")),
            br(),
            downloadButton("btn_download", "Download Sample CSV", class = "btn-block")
          ),
          box(title = "Selected Sample (head)", width = 9, status = "info",
            DTOutput("tbl_sample")
          )
        ),
        fluidRow(
          box(title = "Sample Summary", width = 12,
            verbatimTextOutput("txt_sample_summary")
          )
        )
      ),

      # ── TAB 9: SCRIPT ────────────────────────────────────────────────────────
      tabItem(tabName = "tab_script",
        fluidRow(
          box(title = "Generated R Script", width = 12, status = "info",
            p("This section shows the reproducible R script generated from the workflow executed in the app."),
            verbatimTextOutput("txt_script_preview")
          )
        ),
        fluidRow(
          box(title = "Download Generated Script", width = 4, status = "primary", solidHeader = TRUE,
            downloadButton("btn_download_script", "Download R Script", class = "btn-block")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# ─── SERVER ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  rv <- reactiveValues(
    raw_df      = NULL,
    raw_cens_df = NULL,      # uploaded raw take-all subset
    frame       = NULL,
    framecens   = NULL,      # take-all units (certainty stratum)
    framesamp   = NULL,      # sampling units (frame minus take-all)
    strata      = NULL,
    errors      = NULL,
    solution    = NULL,
    outstrata   = NULL,
    framenew    = NULL,      # sol$framenew — frame with LABEL column for summaryStrata
    x_vars_used = NULL,      # actual X column names passed to buildFrameDF (may include _bin cols)
    eval_res         = NULL,
    sample_df        = NULL,
    optim_plot_files = NULL,  # character vector of per-domain PNG paths
    optim_plot_idx   = 1L,    # currently displayed page
    script_steps     = character(0),
    main_data_expr   = NULL,
    takeall_data_expr = NULL,
    cv_uploads       = NULL   # named list: cv_{i}_dom_{d} -> value, from uploaded CSV
  )

  append_script_step <- function(...) {
    txt <- paste0(...)
    rv$script_steps <- c(rv$script_steps, txt)
  }

  quote_r_string <- function(x) {
    paste0('"', gsub('([\\\\"])', '\\\\\\1', as.character(x)), '"')
  }

  vec_to_r <- function(x) {
    if (length(x) == 0) return('character(0)')
    paste0('c(', paste(vapply(x, quote_r_string, character(1)), collapse = ', '), ')')
  }

  named_numvec_to_r <- function(x) {
    if (length(x) == 0) return('numeric(0)')
    vals <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
    nms  <- names(x)
    if (is.null(nms) || all(is.na(nms)) || all(nms == '')) {
      return(paste0('c(', paste(vals, collapse = ', '), ')'))
    }
    paste0('c(', paste(paste0(quote_r_string(nms), ' = ', vals), collapse = ', '), ')')
  }

  # Bethel baseline computed domain by domain (as the optimiser does) and summed.
  # strata: output of buildStrataDF; errors: data.frame with CVi and domainvalue.
  bethel_total <- function(strata, errors) {
    cv_cols <- grep("^CV", names(errors), value = TRUE)
    doms    <- sort(unique(strata$DOM1))
    sum(vapply(doms, function(d) {
      e   <- errors[as.character(errors$domainvalue) == as.character(d), , drop = FALSE]
      one <- cbind(data.frame(DOM = "DOM1", stringsAsFactors = FALSE),
                   e[1, cv_cols, drop = FALSE])
      sum(ceiling(bethel(strata[strata$DOM1 == d, , drop = FALSE], one)))
    }, numeric(1)))
  }

  # ── HELP MODAL ───────────────────────────────────────────────────────────────
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = tags$span(icon("question-circle"), " SamplingStrata — Workflow Guide"),
      size  = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),

      tabsetPanel(

        # ── Step 1 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("database"), " 1. Data"),
          br(),
          tags$h4("Purpose"),
          tags$p("Load the population dataset that will be used throughout the workflow."),
          tags$h4("Options"),
          tags$ul(
            tags$li(tags$b("Demo dataset:"), " Uses the built-in ",
                    tags$code("swissmunicipalities"), " dataset (2,896 Swiss municipalities
                    with socio-economic variables). Ideal for trying the app without your own data."),
            tags$li(tags$b("Upload CSV:"), " Upload any delimited text file.
                    Choose the separator (comma, semicolon, or tab) and whether the first row
                    is a header.")
          ),
          tags$h4("Action"),
          tags$p("Click ", tags$b("Load Data"), ". A preview (first 100 rows) and a summary
                 will appear on the right.")
        ),

        # ── Step 2 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("sliders-h"), " 2. Variables"),
          br(),
          tags$h4("Purpose"),
          tags$p("Tell SamplingStrata which columns play which role and, optionally,
                 discretise continuous auxiliary variables into bins before building strata."),
          tags$h4("Roles to assign"),
          tags$ul(
            tags$li(tags$b("ID column:"), " Unique record identifier (e.g., municipality code)."),
            tags$li(tags$b("Domain column:"), " Column that defines estimation domains
                    (e.g., region). The optimisation is run separately within each domain."),
            tags$li(tags$b("Target Y variables:"), " The survey variables you want to
                    estimate (e.g., income, population count). Precision constraints will be
                    set for these."),
            tags$li(tags$b("Stratification X variables:"), " Auxiliary variables available
                    in the frame that are correlated with Y and will be used to define strata.")
          ),
          tags$h4("Per-variable categorisation table"),
          tags$p("Once you have ticked one or more X variables, a small table appears below
                 the X selection list with one row per variable:"),
          tags$ul(
            tags$li(tags$b("Categorise? (checkbox)"), " — tick to discretise that variable
                    into equal-frequency bins using ", tags$code("var.bin()"), "."),
            tags$li(tags$b("Bins (numeric)"), " — appears only when the checkbox is ticked;
                    sets the number of bins (2–50). The default is 5. If the variable already
                    has fewer unique values than the requested bins, it is left unchanged and
                    a warning is shown.")
          ),
          tags$p("If an X variable is also selected as a Y target, a new column named ",
                 tags$code("variable_bin"), " is created so that the original column (used
                 as Y) remains untouched."),
          tags$h4("Frequency table"),
          tags$p("After clicking ", tags$b("Prepare Variables"), ", a frequency table is
                 produced for each categorised variable showing the bin label, the count of
                 population units falling in that bin, and the percentage. This helps you
                 check whether the binning produces reasonably balanced groups."),
          tags$h4("Correlation matrix"),
          tags$p("The right-hand panel shows the Pearson correlation between each X and
                 each Y. Strong correlations (close to ±1) indicate that the X variable is
                 a good auxiliary for stratification on the corresponding Y.")
        ),

        # ── Step 3 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("layer-group"), " 3. Frame & Strata"),
          br(),
          tags$h4("Purpose"),
          tags$p("Construct the two main SamplingStrata objects from your variable assignments."),
          tags$h4("Build Frame"),
          tags$p("Calls ", tags$code("buildFrameDF()"), " to create the sampling frame — a
                 standardised data frame with columns ",
                 tags$code("ID"), ", ", tags$code("domainvalue"), ", ",
                 tags$code("Y1…Yk"), ", ", tags$code("X1…Xp"), "."),
          tags$h4("Build Strata"),
          tags$p("Calls ", tags$code("buildStrataDF()"), " to aggregate the frame into
                 strata, computing within-stratum means (", tags$code("M"), ") and
                 standard deviations (", tags$code("S"), ") for each Y variable, plus the
                 stratum population size ", tags$code("N"), "."),
          tags$h4("Check Input"),
          tags$p("Available in Step 4 once the constraints are set: runs ",
                 tags$code("checkInput()"), " to validate the errors, strata and frame
                 objects (no missing values, coherent sizes, etc.) before optimisation."),
          tags$h4("Tip"),
          tags$p("Always run Check Input before proceeding — it surfaces data quality
                 issues early.")
        ),

        # ── Step 4 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("bullseye"), " 4. Constraints"),
          br(),
          tags$h4("Purpose"),
          tags$p("Define the precision targets the optimal allocation must satisfy."),
          tags$h4("CV constraints"),
          tags$p("For each combination of target variable and domain, enter the maximum
                 acceptable ", tags$b("Coefficient of Variation (CV)"), " — expressed as
                 a proportion (e.g., 0.05 = 5 %). The optimiser will search for the
                 smallest sample that keeps all CVs below these thresholds."),
          tags$h4("Bethel Baseline"),
          tags$p("Clicking ", tags$b("Bethel Baseline"), " computes the classical
                 Bethel (multivariate Neyman) allocation on the ", tags$em("current"), "
                 strata, separately for each domain (using that domain's CVs), and sums
                 the results. This gives a lower-bound reference for the optimal sample
                 size before any boundary reshaping."),
          tags$h4("Tip"),
          tags$p("Start with conservative CVs (0.10–0.15). If the resulting sample size
                 is too large, relax the constraints or increase the number of strata
                 in Step 6.")
        ),

        # ── Step 5 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("filter"), " 5. Take-all Strata"),
          br(),
          tags$h4("Purpose"),
          tags$p("This step is ", tags$b("optional"), ". Use it only when some units must be",
                 " included in the sample with certainty (for example very large firms,",
                 " hospitals, or municipalities). These units form the take-all or census stratum."),
          tags$h4("How to use"),
          tags$ol(
            tags$li(tags$b("Enable take-all strata"), " — tick the checkbox only if you want",
                    " to use an external list of certainty units."),
            tags$li(tags$b("Upload a CSV file"), " — the uploaded file must be a subset of the",
                    " dataset loaded in Step 1."),
            tags$li(tags$b("Validation"), " — the app checks that the uploaded file contains only",
                    " units that belong to the original data and that the required variables are present."),
            tags$li(tags$b("Build take-all frame"), " — if the file is valid, the same ",
                    tags$code("buildFrameDF()"), " settings chosen in Step 2 are applied to the uploaded subset.",
                    " The resulting object becomes ", tags$code("framecens"), "."),
            tags$li(tags$b("Sampling frame"), " — the remaining units from the full frame become ",
                    tags$code("framesamp"), ".")
          ),
          tags$h4("What changes downstream"),
          tags$p("When take-all strata are active, the optimisation passes both ",
                 tags$code("framesamp"), " and ", tags$code("framecens"),
                 " to ", tags$code("optimStrata()"), ". The take-all units contribute",
                 " to the precision of estimates but are not subject to boundary reshaping."),
          tags$h4("Tip"),
          tags$p("If you rebuild the frame (Step 3) or change the variable settings (Step 2),",
                 " upload and validate the take-all subset again.")
        ),

        # ── Step 6 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("cogs"), " 6. Optimization"),
          br(),
          tags$h4("Purpose"),
          tags$p("Find the optimal strata and sample allocation that
                 minimise total sample size while meeting all CV constraints."),
          tags$h4("Method"),
          tags$ul(
            tags$li(tags$b("atomic:"), " Atomic method
                    (", tags$code("method='atomic'"), ") — works with 
                    both categorical or continuous stratification variables."),
            tags$li(tags$b("continuous:"), " Continuous method
                    (", tags$code("method='continuous'"), ") — faster but
                    works only with continuous auxiliary variables.")
          ),
          tags$h4("Key parameters"),
          tags$ul(
            tags$li(tags$b("Max strata per domain:"), " Upper bound on the number of
                    strata the algorithm may create per domain."),
            tags$li(tags$b("Use nStrata from KmeansSolution:"), " Strongly recommended
                    to choose the number of aggregated strata in each domain"),            
            tags$li(tags$b("Iterations (iter):"), " Number of generations for the
                    genetic algorithm. More iterations → better solution, longer runtime."),
            tags$li(tags$b("Population size (pops):"), " Number of candidate solutions
                    evolved in parallel."),
            tags$li(tags$b("Use KmeansSolution as seed:"), " Pre-initialises the GA
                    with a K-means clustering — strongly recommended, speeds up
                    convergence - to be used with atomic method."),
            tags$li(tags$b("Use suggestions (continuous method only):"), " Runs ",
                    tags$code("suggestions <- prepareSuggestion(frame, nStrata)"),
                    " and passes ", tags$code("suggestions"), " to ",
                    tags$code("optimStrata()"), " as the starting solution.")
          ),
          
          tags$h4("Output"),
          tags$p("The convergence plot shows how the objective (sample size) decreases
                 across iterations. The log and strata summary table confirm the final
                 allocation and expected CVs.")
        ),

        # ── Step 7 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("chart-bar"), " 7. Results"),
          br(),
          tags$h4("Purpose"),
          tags$p("Evaluate and, if needed, rescale the optimised solution."),
          tags$h4("Evaluate Solution"),
          tags$p("Calls ", tags$code("evalSolution()"), " which draws ",
                 tags$b("nsampl"), " repeated samples from the frame and computes
                 empirical CVs and bias for each target variable and domain.
                 The resulting plots show whether the solution really meets the
                 precision constraints in practice."),
          tags$h4("Adjust Sample Size"),
          tags$p("If the optimised total differs from an operational target (e.g., a
                 budget constraint fixes n = 250), enter the desired total and click
                 ", tags$b("Adjust Size"), ". The per-stratum allocations are rescaled
                 proportionally (minimum 2 units per stratum; each stratum is also
                 capped at its population size N)."),
          tags$h4("Strata Visualisation"),
          tags$p("Use the X/Y axis dropdowns and the Domain selector to inspect
                 where the strata boundaries fall in the auxiliary-variable space.")
        ),

        # ── Step 8 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("table"), " 8. Sample"),
          br(),
          tags$h4("Purpose"),
          tags$p("Draw the actual sample units from the frame according to the
                 optimised (or adjusted) allocation."),
          tags$h4("Method"),
          tags$ul(
            tags$li(tags$b("Random (selectSample):"), " Simple random sampling
                    without replacement within each stratum — the standard choice."),
            tags$li(tags$b("Systematic (selectSampleSystematic):"), " Systematic
                    sampling within strata. Requires a sort variable; useful when
                    implicit stratification along an ordering variable is desired.")
          ),
          tags$h4("Sort variable"),
          tags$p("Only relevant for the systematic method. Choose the variable by
                 which units are ordered within each stratum before systematic
                 selection."),
          tags$h4("Download"),
          tags$p("After selection, click ", tags$b("Download Sample CSV"), " to save
                 the selected units with all their frame variables plus the stratum
                 label.")
        ),

        # ── Step 9 ──────────────────────────────────────────────────────────
        tabPanel(tags$span(icon("file-code"), " 9. R Script"),
          br(),
          tags$h4("Purpose"),
          tags$p("Review and download a reproducible R script generated from the
                 workflow executed in the app."),
          tags$h4("Preview"),
          tags$p("The full script is shown on screen before download so you can
                 inspect the recorded commands."),
          tags$h4("Content"),
          tags$p("The script includes an initial description of what it does, the
                 expected input files, and the output files it may create. It then
                 reproduces the executable SamplingStrata workflow recorded during
                 the interaction."),
          tags$h4("Download"),
          tags$p("Click ", tags$b("Download R Script"), " to save the generated
                 script as a standalone .R file.")
        )

      ) # end tabsetPanel
    )) # end modalDialog / showModal
  })

  # ── STEP 1: Load data ───────────────────────────────────────────────────────
  observeEvent(input$btn_load, {
    if (input$data_source == "demo") {
      data("swissmunicipalities", package = "SamplingStrata")
      rv$raw_df <- swissmunicipalities
      rv$main_data_expr <- 'data("swissmunicipalities", package = "SamplingStrata")
raw_df <- swissmunicipalities'
    } else {
      req(input$upload_file)
      rv$raw_df <- read.csv(input$upload_file$datapath,
                            header = input$header,
                            sep    = input$sep)
      rv$main_data_expr <- paste0(
        'raw_df <- read.csv(', quote_r_string(input$upload_file$name),
        ', header = ', if (isTRUE(input$header)) 'TRUE' else 'FALSE',
        ', sep = ', quote_r_string(input$sep), ')'
      )
    }
    # Reset everything downstream of the raw data to avoid stale objects
    for (nm in c("raw_cens_df", "frame", "framecens", "framesamp", "strata",
                 "errors", "solution", "outstrata", "framenew", "x_vars_used",
                 "eval_res", "sample_df", "optim_plot_files",
                 "takeall_data_expr", "cv_uploads")) {
      rv[[nm]] <- NULL
    }
    rv$optim_plot_idx <- 1L
    frame_xbins(NULL)
    xfreq_df(NULL)
    check_msg("Set constraints first (Tab 4), then click Check Input.")
    rv$script_steps <- c(
      '# Reproducible script generated from the SamplingStrata Shiny app',
      'library(SamplingStrata)',
      '',
      '# 1. Load main data',
      rv$main_data_expr,
      ''
    )
    showNotification("Data loaded.", type = "message")
  })

  output$tbl_data <- renderDT({
    req(rv$raw_df)
    datatable(head(rv$raw_df, 100), options = list(scrollX = TRUE, pageLength = 8))
  })

  output$txt_summary <- renderPrint({
    req(rv$raw_df)
    summary(rv$raw_df)
  })

  # Message shown in the "Check Input Result" box (updated by btn_check)
  check_msg <- reactiveVal("Set constraints first (Tab 4), then click Check Input.")

  # ── STEP 2: Variable assignment ─────────────────────────────────────────────
  col_choices <- reactive({
    req(rv$raw_df)
    names(rv$raw_df)
  })

  # Helper: case-insensitive match of `wanted` names against the actual columns
  # of the loaded dataset. Returns the actual column names (preserving the
  # original capitalisation) for those that exist, dropping the rest.
  match_cols_ci <- function(wanted, cols) {
    if (is.null(cols) || length(cols) == 0) return(character(0))
    idx <- match(tolower(wanted), tolower(cols))
    cols[idx[!is.na(idx)]]
  }

  # Detect whether the currently loaded dataset looks like the swissmunicipalities
  # demo (case-insensitive on column names), so we can pre-select sensible
  # defaults in the variable assignment UI.
  is_swiss_demo <- reactive({
    req(rv$raw_df)
    if (!isTRUE(input$data_source == "demo")) return(FALSE)
    expected <- c("id", "Reg", "Airbat", "Airind",
                  "Surfacesbois", "Surfacescult", "HApoly", "POPTOT")
    length(match_cols_ci(expected, names(rv$raw_df))) == length(expected)
  })

  output$ui_id_var <- renderUI({
    sel <- if (isTRUE(is_swiss_demo()))
      match_cols_ci("id", col_choices()) else NULL
    selectInput("id_var", "ID column",
                choices = col_choices(), selected = sel)
  })
  output$ui_domain_var <- renderUI({
    sel <- if (isTRUE(is_swiss_demo()))
      match_cols_ci("Reg", col_choices()) else NULL
    selectInput("domain_var", "Domain column",
                choices = col_choices(), selected = sel)
  })
  output$ui_y_vars <- renderUI({
    sel <- if (isTRUE(is_swiss_demo())) match_cols_ci(
      c("Airbat", "Airind", "Surfacesbois", "Surfacescult"),
      col_choices()) else NULL
    checkboxGroupInput("y_vars", "Target (Y) variables",
                       choices = col_choices(), selected = sel, inline = FALSE)
  })
  output$ui_x_vars <- renderUI({
    sel <- if (isTRUE(is_swiss_demo()))
      match_cols_ci(c("HApoly", "POPTOT"), col_choices()) else NULL
    checkboxGroupInput("x_vars", "Stratification (X) variables",
                       choices = col_choices(), selected = sel, inline = FALSE)
  })

  # Per-variable categorisation table (shown once X vars are selected)
  output$ui_x_bins_table <- renderUI({
    vars <- input$x_vars
    if (is.null(vars) || length(vars) == 0) return(NULL)
    tagList(
      tags$hr(),
      tags$p(tags$b("Categorisation settings (per X variable):")),
      tags$table(
        class = "table table-condensed table-bordered",
        style = "font-size:12px; margin-bottom:6px;",
        tags$thead(tags$tr(
          tags$th("Variable"),
          tags$th(style = "text-align:center;", "Categorise?"),
          tags$th("Bins")
        )),
        tags$tbody(
          lapply(vars, function(v) {
            tags$tr(
              tags$td(v),
              tags$td(style = "text-align:center;",
                checkboxInput(paste0("cat_", v), label = NULL, value = FALSE)
              ),
              tags$td(
                conditionalPanel(
                  condition = sprintf("input['cat_%s'] == true", v),
                  numericInput(paste0("bins_", v), label = NULL,
                               value = 5, min = 2, max = 50, width = "80px")
                )
              )
            )
          })
        )
      )
    )
  })

  output$plot_corr <- renderPlot({
    req(rv$raw_df, input$x_vars, input$y_vars)
    xv <- input$x_vars
    yv <- input$y_vars
    df <- rv$raw_df
    # build numeric correlation matrix X vs Y
    tryCatch({
      mat <- cor(df[, xv, drop = FALSE], df[, yv, drop = FALSE], use = "complete.obs")
      corrplot::corrplot(mat, method = "color", addCoef.col = "black",
                         tl.col = "black", cl.ratio = 0.2,
                         title = "Correlation X vs Y", mar = c(0,0,2,0))
    }, error = function(e) {
      plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 1.2)
    })
  })

  # Apply per-variable var.bin
  frame_xbins <- reactiveVal(NULL)
  xfreq_df    <- reactiveVal(NULL)

  observeEvent(input$btn_prepare, {
    req(rv$raw_df, input$x_vars, input$y_vars, input$id_var, input$domain_var)
    df_orig   <- rv$raw_df
    df        <- df_orig
    x_vars    <- input$x_vars
    any_cat   <- FALSE
    freq_list <- list()

    withProgress(message = "Preparing variables...", {
      for (i in seq_along(x_vars)) {
        v      <- x_vars[i]
        do_cat <- isTRUE(input[[paste0("cat_", v)]])
        nb     <- if (do_cat) input[[paste0("bins_", v)]] else 0L

        if (do_cat && !is.null(nb) && nb >= 2) {
          col <- df_orig[[v]]
          if (length(unique(col)) > nb) {
            any_cat <- TRUE
            if (v %in% input$y_vars) {
              # X is also a Y: create a _bin copy to leave Y column untouched
              bin_name       <- paste0(v, "_bin")
              df[[bin_name]] <- var.bin(col, nb)
              x_vars[i]      <- bin_name
              ft             <- as.data.frame(table(Category = df[[bin_name]]))
              ft$Variable    <- paste0(v, " \u2192 ", bin_name)
            } else {
              df[[v]]     <- var.bin(col, nb)
              ft          <- as.data.frame(table(Category = df[[v]]))
              ft$Variable <- v
            }
            ft$Percent <- round(100 * ft$Freq / sum(ft$Freq), 1)
            freq_list[[length(freq_list) + 1]] <-
              ft[, c("Variable", "Category", "Freq", "Percent")]
          } else {
            showNotification(
              paste0(v, ": fewer unique values than requested bins — skipped."),
              type = "warning", duration = 6)
          }
        }
      }
    })

    if (any_cat) {
      showNotification("Selected X variables categorised.", type = "message")
      xfreq_df(do.call(rbind, freq_list))
    } else {
      showNotification("No categorisation applied.", type = "message")
      xfreq_df(NULL)
    }

    frame_xbins(df)
    rv$x_vars_used <- x_vars    # store actual X names (may include _bin columns)

    prep_lines <- c(
      '# 2. Prepare variables',
      'df_prepared <- raw_df'
    )
    for (v in input$x_vars) {
      do_cat <- isTRUE(input[[paste0("cat_", v)]])
      nb     <- if (do_cat) input[[paste0("bins_", v)]] else NULL
      if (isTRUE(do_cat) && !is.null(nb) && nb >= 2) {
        if (v %in% input$y_vars) {
          prep_lines <- c(prep_lines,
            paste0('df_prepared[[', quote_r_string(paste0(v, '_bin')), ']] <- var.bin(df_prepared[[', quote_r_string(v), ']], ', nb, ')'))
        } else {
          prep_lines <- c(prep_lines,
            paste0('df_prepared[[', quote_r_string(v), ']] <- var.bin(df_prepared[[', quote_r_string(v), ']], ', nb, ')'))
        }
      }
    }
    if (!length(input$x_vars)) prep_lines <- c(prep_lines, '# No X variables selected')
    append_script_step(paste(prep_lines, collapse = '
'), '
')
  })

  output$tbl_xbins <- renderDT({
    req(frame_xbins(), rv$x_vars_used)
    show_cols <- intersect(c(input$id_var, rv$x_vars_used), names(frame_xbins()))
    datatable(head(frame_xbins()[, show_cols, drop = FALSE], 50),
              options = list(scrollX = TRUE, pageLength = 8))
  })

  output$ui_xfreq_note <- renderUI({
    if (is.null(xfreq_df()))
      tags$p(tags$em(
        "No variables were categorised — tick the 'Categorise?' box for one or more X variables
         and click Prepare Variables to see the frequency table here."
      ))
    else NULL
  })

  output$tbl_xfreq <- renderDT({
    req(xfreq_df())
    datatable(
      xfreq_df(),
      colnames  = c("Variable", "Category", "Count", "%"),
      rownames  = FALSE,
      options   = list(scrollX = TRUE, pageLength = 20)
    )
  })

  # ── STEP 3: Frame & Strata ──────────────────────────────────────────────────
  observeEvent(input$btn_build_frame, {
    req(input$id_var, input$x_vars, input$y_vars, input$domain_var)
    if (is.null(frame_xbins())) {
      showNotification("Click 'Prepare Variables' (Step 2) before building the frame.",
                       type = "warning", duration = 8)
      return()
    }
    x_vars_final <- if (!is.null(rv$x_vars_used)) rv$x_vars_used else input$x_vars
    tryCatch({
      rv$frame <- buildFrameDF(
        df          = frame_xbins(),
        id          = input$id_var,
        X           = x_vars_final,      # actual X names (may include _bin cols)
        Y           = input$y_vars,
        domainvalue = input$domain_var
      )
      append_script_step(
        '# 3. Build frame
',
        'frame <- buildFrameDF(df = df_prepared, id = ', quote_r_string(input$id_var),
        ', X = ', vec_to_r(x_vars_final),
        ', Y = ', vec_to_r(input$y_vars),
        ', domainvalue = ', quote_r_string(input$domain_var), ')
'
      )
      showNotification("Frame built.", type = "message")
    }, error = function(e) {
      showNotification(paste("Frame error:", e$message), type = "error", duration = 10)
    })
  })

  output$tbl_frame <- renderDT({
    req(rv$frame)
    datatable(head(rv$frame, 100), options = list(scrollX = TRUE, pageLength = 8))
  })

  observeEvent(input$btn_build_strata, {
    req(rv$frame)
    tryCatch({
      rv$strata <- buildStrataDF(rv$frame)
      append_script_step('# 4. Build strata
strata <- buildStrataDF(frame)
')
      showNotification("Strata built.", type = "message")
    }, error = function(e) {
      showNotification(paste("Strata error:", e$message), type = "error", duration = 10)
    })
  })

  output$tbl_strata <- renderDT({
    req(rv$strata)
    datatable(head(rv$strata, 100), options = list(scrollX = TRUE, pageLength = 8))
  })

  output$txt_check <- renderPrint({ cat(check_msg(), "\n") })

  observeEvent(input$btn_check, {
    if (is.null(rv$strata) || is.null(rv$frame) || is.null(rv$errors)) {
      check_msg("Build frame and strata (Tab 3) and set constraints (Tab 4) first.")
      return()
    }
    tryCatch({
      checkInput(errors = rv$errors, strata = rv$strata, sampframe = rv$frame)
      check_msg("checkInput: OK")
      append_script_step('# Optional: check input
checkInput(errors = errors, strata = strata, sampframe = frame)
')
      showNotification("checkInput passed.", type = "message")
    }, error = function(e) {
      check_msg(paste("ERROR:", e$message))
      showNotification(paste("Check failed:", e$message), type = "error", duration = 10)
    }, warning = function(w) {
      check_msg(paste("WARNING:", w$message))
      showNotification(paste("Check warning:", w$message), type = "warning", duration = 10)
    })
  })

  # ── STEP 4: Constraints ──────────────────────────────────────────────────────

  # Global default CV inputs (one per Y variable, left panel)
  output$ui_cv_global <- renderUI({
    req(input$y_vars)
    yv <- input$y_vars
    tagList(
      lapply(seq_along(yv), function(i) {
        numericInput(
          inputId = paste0("cv_global_", i),
          label   = paste("Default CV for", yv[i]),
          value   = 0.10, min = 0.001, max = 1, step = 0.01
        )
      })
    )
  })

  # "Apply to all domains" — push global defaults into every per-domain input
  observeEvent(input$btn_apply_global_cv, {
    req(input$y_vars, rv$strata)
    yv   <- input$y_vars
    doms <- sort(unique(rv$strata$DOM1))
    for (i in seq_along(yv)) {
      global_val <- input[[paste0("cv_global_", i)]]
      if (is.null(global_val) || is.na(global_val)) global_val <- 0.10
      for (d in seq_along(doms)) {
        updateNumericInput(session, paste0("cv_", i, "_dom_", d), value = global_val)
      }
    }
    showNotification("Global CV applied to all domains.", type = "message")
  })

  # Per-domain CV inputs (right panel)
  output$ui_errors_form <- renderUI({
    req(input$y_vars, rv$strata)
    yv      <- input$y_vars
    doms    <- sort(unique(rv$strata$DOM1))
    uploads <- rv$cv_uploads   # may be NULL

    tagList(
      lapply(seq_along(yv), function(i) {
        box_content <- lapply(seq_along(doms), function(d) {
          key <- paste0("cv_", i, "_dom_", d)
          val <- if (!is.null(uploads) && !is.null(uploads[[key]])) uploads[[key]] else 0.10
          numericInput(
            inputId = key,
            label   = paste0("Dom ", doms[d]),
            value   = val, min = 0.001, max = 1, step = 0.01
          )
        })
        tagList(
          tags$hr(),
          tags$b(paste("Max CV for", yv[i])),
          tags$div(
            style = "display:flex; flex-wrap:wrap; gap:6px;",
            box_content
          )
        )
      })
    )
  })

  observeEvent(input$btn_set_errors, {
    req(rv$strata, input$y_vars)
    yv   <- input$y_vars
    doms <- sort(unique(rv$strata$DOM1))

    # Build errors data.frame: one row per domain, one CVi column per Y variable
    errs <- data.frame(DOM = rep("DOM1", length(doms)))
    for (i in seq_along(yv)) {
      errs[[paste0("CV", i)]] <- sapply(seq_along(doms), function(d) {
        v <- input[[paste0("cv_", i, "_dom_", d)]]
        if (is.null(v) || is.na(v)) 0.10 else v
      })
    }
    errs$domainvalue <- doms   # domainvalue must be the LAST column
    rv$errors <- errs

    # Build script lines reflecting per-domain CVs
    cv_lines <- vapply(seq_along(yv), function(i) {
      cv_vec <- sapply(seq_along(doms), function(d) {
        v <- input[[paste0("cv_", i, "_dom_", d)]]
        if (is.null(v) || is.na(v)) 0.10 else v
      })
      paste0('errors$CV', i, ' <- c(', paste(format(cv_vec, scientific = FALSE), collapse = ', '), ')')
    }, character(1))

    append_script_step(
      '# 5. Set precision constraints
',
      'errors <- data.frame(DOM = rep("DOM1", ', length(doms), '))
',
      paste(cv_lines, collapse = '\n'),
      '\nerrors$domainvalue <- ',
      if (is.numeric(doms)) paste0('c(', paste(doms, collapse = ', '), ')')
      else vec_to_r(as.character(doms)),
      '\n'
    )
    showNotification("Constraints set.", type = "message")
  })

  # Bethel baseline: see bethel_total() (computed per domain, then summed).

  output$tbl_errors <- renderDT({
    req(rv$errors)
    datatable(rv$errors, options = list(scrollX = TRUE))
  })

  output$txt_bethel <- renderPrint({
    req(rv$strata, rv$errors)
    tryCatch({
      cat("Bethel total sample size (sum over domains):",
          bethel_total(rv$strata, rv$errors), "\n")
    }, error = function(e) cat("Error:", e$message, "\n"))
  })

  observeEvent(input$btn_bethel, {
    req(rv$strata, rv$errors)
    tryCatch({
      n_bethel <- bethel_total(rv$strata, rv$errors)
      append_script_step('# Optional: Bethel baseline (domain by domain, then summed)
cv_cols <- grep("^CV", names(errors), value = TRUE)
n_bethel <- sum(sapply(sort(unique(strata$DOM1)), function(d) {
  e <- errors[as.character(errors$domainvalue) == as.character(d), , drop = FALSE]
  one <- cbind(data.frame(DOM = "DOM1"), e[1, cv_cols, drop = FALSE])
  sum(ceiling(bethel(strata[strata$DOM1 == d, , drop = FALSE], one)))
}))
n_bethel
')
      showNotification(paste("Bethel n =", n_bethel),
                       type = "message", duration = 8)
    }, error = function(e) showNotification(paste("Bethel error:", e$message),
                                             type = "error", duration = 10))
  })

  # ── CV download ──────────────────────────────────────────────────────────────
  # Build a tidy CSV: variable, domainvalue, cv  (one row per variable × domain)
  output$btn_download_cv <- downloadHandler(
    filename = function() paste0("cv_constraints_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      req(input$y_vars, rv$strata)
      yv   <- input$y_vars
      doms <- sort(unique(rv$strata$DOM1))
      rows <- do.call(rbind, lapply(seq_along(yv), function(i) {
        do.call(rbind, lapply(seq_along(doms), function(d) {
          v <- input[[paste0("cv_", i, "_dom_", d)]]
          data.frame(variable    = yv[i],
                     domainvalue = doms[d],
                     cv          = if (is.null(v) || is.na(v)) 0.10 else v,
                     stringsAsFactors = FALSE)
        }))
      }))
      write.csv(rows, file, row.names = FALSE)
    }
  )

  # ── CV upload ────────────────────────────────────────────────────────────────
  output$ui_cv_upload_status <- renderUI({ NULL })

  observeEvent(input$upload_cv, {
    req(input$upload_cv, input$y_vars, rv$strata)
    yv   <- input$y_vars
    doms <- sort(unique(rv$strata$DOM1))

    df <- tryCatch(
      read.csv(input$upload_cv$datapath, stringsAsFactors = FALSE),
      error = function(e) e
    )
    if (inherits(df, "error")) {
      output$ui_cv_upload_status <- renderUI(
        tags$p(style = "color:red;", icon("times-circle"), " Could not read CSV: ", df$message)
      )
      return()
    }
    required <- c("variable", "domainvalue", "cv")
    if (!all(required %in% names(df))) {
      output$ui_cv_upload_status <- renderUI(
        tags$p(style = "color:red;", icon("times-circle"),
               " CSV must have columns: variable, domainvalue, cv")
      )
      return()
    }

    # Build named list: cv_{i}_dom_{d} -> value
    new_vals <- list()
    warnings <- character(0)
    for (i in seq_along(yv)) {
      for (d in seq_along(doms)) {
        row <- df[df$variable == yv[i] & as.character(df$domainvalue) == as.character(doms[d]), ]
        if (nrow(row) == 0) {
          warnings <- c(warnings, paste0(yv[i], " / dom ", doms[d]))
          new_vals[[paste0("cv_", i, "_dom_", d)]] <- 0.10
        } else {
          new_vals[[paste0("cv_", i, "_dom_", d)]] <- as.numeric(row$cv[1])
        }
      }
    }
    rv$cv_uploads <- new_vals   # triggers ui_errors_form to redraw with new values

    if (length(warnings) > 0) {
      output$ui_cv_upload_status <- renderUI(
        tags$p(style = "color:orange;", icon("exclamation-triangle"),
               paste(" Defaulted to 0.10 for:", paste(warnings, collapse = "; ")))
      )
    } else {
      output$ui_cv_upload_status <- renderUI(
        tags$p(style = "color:green;", icon("check-circle"),
               " CV values loaded successfully.")
      )
    }
    showNotification("CV constraints uploaded.", type = "message")
  })

  # ── STEP 5: Take-all strata ──────────────────────────────────────────────────

  # Reset take-all objects whenever the main frame is rebuilt
  observeEvent(rv$frame, {
    rv$raw_cens_df <- NULL
    rv$framecens   <- NULL
    rv$framesamp   <- NULL
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$use_cens, {
    if (!isTRUE(input$use_cens)) {
      rv$raw_cens_df <- NULL
      rv$framecens   <- NULL
      rv$framesamp   <- NULL
      rv$takeall_data_expr <- NULL
    }
  }, ignoreInit = TRUE)

  observeEvent(input$btn_apply_cens, {
    req(isTRUE(input$use_cens), rv$raw_df, rv$frame, input$cens_file,
        input$id_var, input$domain_var, input$y_vars)

    subset_df <- tryCatch(
      read.csv(input$cens_file$datapath,
               header = input$cens_header,
               sep    = input$cens_sep,
               stringsAsFactors = FALSE),
      error = function(e) e
    )

    if (inherits(subset_df, "error")) {
      showNotification(paste("Take-all file error:", subset_df$message),
                       type = "error", duration = 10)
      return()
    }

    x_vars_final <- if (!is.null(rv$x_vars_used)) rv$x_vars_used else input$x_vars
    req(length(x_vars_final) > 0)

    if (is.null(frame_xbins())) {
      showNotification("Click 'Prepare Variables' (Step 2) first.", type = "warning")
      return()
    }

    # Only the original (un-binned) columns must be in the uploaded file: the
    # binned versions are taken from the prepared data, so that framecens uses
    # exactly the same categorisation as the main frame.
    required_cols <- unique(c(input$id_var, input$domain_var, input$y_vars, input$x_vars))
    missing_cols  <- setdiff(required_cols, names(subset_df))
    if (length(missing_cols) > 0) {
      showNotification(
        paste0("Invalid take-all file: missing required column(s): ",
               paste(missing_cols, collapse = ", ")),
        type = "error", duration = 12)
      return()
    }

    raw_ids <- as.character(rv$raw_df[[input$id_var]])
    if (anyDuplicated(raw_ids)) {
      showNotification(
        paste0("The ID column '", input$id_var,
               "' in the main data is not unique. Subset validation requires unique IDs."),
        type = "error", duration = 12)
      return()
    }

    subset_ids <- as.character(subset_df[[input$id_var]])
    if (length(subset_ids) == 0) {
      showNotification("The uploaded take-all file is empty.", type = "error", duration = 10)
      return()
    }
    if (anyDuplicated(subset_ids)) {
      showNotification(
        paste0("Invalid take-all file: duplicated ID values found in column '",
               input$id_var, "'."),
        type = "error", duration = 12)
      return()
    }

    extra_ids <- setdiff(subset_ids, raw_ids)
    if (length(extra_ids) > 0) {
      showNotification(
        paste0("Invalid take-all file: ", length(extra_ids),
               " ID(s) are not present in the main data. Example: ",
               paste(utils::head(extra_ids, 5), collapse = ", ")),
        type = "error", duration = 15)
      return()
    }

    raw_match_idx <- match(subset_ids, raw_ids)
    compare_cols  <- intersect(required_cols, names(rv$raw_df))
    mismatches <- vapply(compare_cols, function(colname) {
      lhs <- as.character(subset_df[[colname]])
      rhs <- as.character(rv$raw_df[[colname]][raw_match_idx])
      isTRUE(any(is.na(raw_match_idx) | (lhs != rhs) | (is.na(lhs) != is.na(rhs)),
                 na.rm = TRUE))
    }, logical(1))

    if (any(mismatches)) {
      bad_cols <- names(mismatches)[mismatches]
      showNotification(
        paste0("Invalid take-all file: the uploaded rows do not match the main dataset for column(s): ",
               paste(bad_cols, collapse = ", "), "."),
        type = "error", duration = 15)
      return()
    }

    # Take the prepared (binned) rows matching the uploaded IDs
    prep_df     <- frame_xbins()
    cens_prep   <- prep_df[match(subset_ids, as.character(prep_df[[input$id_var]])), ,
                           drop = FALSE]

    built_cens <- tryCatch({
      buildFrameDF(
        df          = cens_prep,
        id          = input$id_var,
        X           = x_vars_final,
        Y           = input$y_vars,
        domainvalue = input$domain_var
      )
    }, error = function(e) e)

    if (inherits(built_cens, "error")) {
      showNotification(paste("Take-all buildFrameDF error:", built_cens$message),
                       type = "error", duration = 12)
      return()
    }

    if (nrow(built_cens) == 0) {
      showNotification("The uploaded take-all file produced an empty framecens.",
                       type = "error", duration = 10)
      return()
    }
    if (nrow(built_cens) >= nrow(rv$frame)) {
      showNotification("All units are in the take-all file: sampling frame would be empty.",
                       type = "error", duration = 10)
      return()
    }

    rv$raw_cens_df <- subset_df
    rv$framecens   <- built_cens
    rv$framesamp   <- rv$frame[!rv$frame$id %in% rv$framecens$id, , drop = FALSE]
    rv$takeall_data_expr <- paste0(
      'takeall_raw <- read.csv(', quote_r_string(input$cens_file$name),
      ', header = ', if (isTRUE(input$cens_header)) 'TRUE' else 'FALSE',
      ', sep = ', quote_r_string(input$cens_sep), ')'
    )
    append_script_step(
      '# 6. Optional take-all subset
',
      rv$takeall_data_expr, '
',
      'takeall_prepared <- df_prepared[match(as.character(takeall_raw[[', quote_r_string(input$id_var),
      ']]), as.character(df_prepared[[', quote_r_string(input$id_var), ']])), , drop = FALSE]
',
      'framecens <- buildFrameDF(df = takeall_prepared, id = ', quote_r_string(input$id_var),
      ', X = ', vec_to_r(x_vars_final),
      ', Y = ', vec_to_r(input$y_vars),
      ', domainvalue = ', quote_r_string(input$domain_var), ')
',
      'framesamp <- frame[!frame$id %in% framecens$id, , drop = FALSE]
'
    )

    all_doms  <- sort(unique(rv$frame$domainvalue))
    cens_doms <- sort(unique(rv$framecens$domainvalue))
    missing_d <- setdiff(all_doms, cens_doms)

    if (length(missing_d) > 0) {
      showNotification(
        paste0("Take-all subset validated: ", nrow(rv$framecens), " units in framecens, ",
               nrow(rv$framesamp), " units in framesamp. ",
               "Some domains have no take-all units (",
               paste(missing_d, collapse = ", "), "), which is allowed."),
        type = "message", duration = 12)
    } else {
      showNotification(
        paste0("Take-all subset validated: ", nrow(rv$framecens), " units in framecens, ",
               nrow(rv$framesamp), " units in framesamp. All ",
               length(all_doms), " domains contain take-all units."),
        type = "message", duration = 12)
    }
  })

  output$txt_cens_preview <- renderPrint({
    if (!isTRUE(input$use_cens)) {
      cat("Take-all strata disabled (optional step).
")
    } else if (is.null(rv$framecens)) {
      cat("Upload a CSV subset and click 'Validate and build take-all frame'.
")
      cat("The uploaded file must be a subset of the main dataset loaded in Step 1.
")
    } else {
      cat("Uploaded take-all raw rows:", nrow(rv$raw_cens_df), "
")
      cat("Take-all units  (framecens):", nrow(rv$framecens), "
")
      cat("Sampling units  (framesamp):", nrow(rv$framesamp), "
")
      cat("Total frame:               ", nrow(rv$frame),     "
")
      cat("Fraction censused:         ",
          round(100 * nrow(rv$framecens) / nrow(rv$frame), 1), "%

")
      all_doms  <- sort(unique(rv$frame$domainvalue))
      cens_doms <- sort(unique(rv$framecens$domainvalue))
      missing_d <- setdiff(all_doms, cens_doms)
      if (length(missing_d) == 0) {
        cat("Domain coverage: OK — all", length(all_doms), "domains have take-all units.
")
      } else {
        cat("Domain coverage: PARTIAL
")
        cat("  Domains with take-all units   :", paste(cens_doms, collapse = ", "), "
")
        cat("  Domains WITHOUT take-all units:", paste(missing_d, collapse = ", "), "
")
        cat("  This is allowed: optimisation can still run.
")
      }
    }
  })

  output$tbl_cens_head <- renderDT({
    req(isTRUE(input$use_cens), rv$framecens)
    datatable(head(rv$framecens, 50),
              options = list(scrollX = TRUE, pageLength = 8),
              caption = "First 50 take-all units")
  })

# ── STEP 6: Optimisation ─────────────────────────────────────────────────────
  # When any X variable has been categorised (xfreq_df is non-NULL after
  # Prepare Variables) only "atomic" is valid; otherwise both are available.
  output$ui_optim_method <- renderUI({
    categorised <- !is.null(xfreq_df())
    if (categorised) {
      tagList(
        selectInput("optim_method", "Method", choices = c("atomic"), selected = "atomic"),
        helpText(icon("info-circle"),
                 " Only 'atomic' is available when one or more X variables are categorised.")
      )
    } else {
      selectInput("optim_method", "Method",
                  choices = c("atomic", "continuous"), selected = "continuous")
    }
  })
  optim_log <- reactiveVal("")

  # Sort convergence plot files by their numeric index (optim_plot_<n>_<time>.png);
  # plain sort() breaks once the index exceeds 999 (e.g. "1000" < "999").
  sort_plot_files <- function(f) {
    if (length(f) == 0) return(f)
    f[order(as.integer(sub("^optim_plot_([0-9]+)_.*$", "\\1", basename(f))))]
  }

  # ── KMeans nStrata (runs when checkbox is activated or strata change) ────────
  kmeans_nstrata_vec    <- reactiveVal(NULL)   # per-domain nStrata vector
  kmeans_solution_stored <- reactiveVal(NULL)  # full KmeansSolution result (for seed)

  # Re-run when checkbox toggled OR when method changes (continuous vs atomic
  # require different KMeans functions: KmeansSolution2 vs KmeansSolution).
  # Signature of the inputs the stored KMeans result was computed from. If strata,
  # errors, method or maxclusters change afterwards, the stored seed is outdated
  # and passing it to optimStrata() fails ("number of items to replace ...").
  kmeans_sig <- reactiveVal(NULL)
  km_sig <- function(method, maxcl) list(method, maxcl, rv$errors, rv$strata)

  run_kmeans <- function() {
    if (!isTRUE(input$kmeans_nstrata)) {
      kmeans_nstrata_vec(NULL)
      kmeans_solution_stored(NULL)
      return()
    }
    method <- if (is.null(input$optim_method)) "atomic" else input$optim_method
    maxcl  <- max(2L, as.integer(if (is.null(input$nstrata)) 10L else input$nstrata))

    if (method == "continuous") {
      req(rv$frame, rv$errors)
      withProgress(message = "Running KmeansSolution2 for nStrata…", {
        tryCatch({
          km  <- KmeansSolution2(frame = rv$frame, errors = rv$errors,
                                 maxclusters = maxcl)
          tbl <- tapply(km$suggestions, km$domainvalue,
                        function(x) length(unique(x)))
          # Convert to a plain named integer vector — tapply returns a 1D array
          # which can cause colnames<- failures deep in optimizeStrata2.
          # setNames on as.integer() drops the dim attribute while keeping names.
          vec <- setNames(as.integer(pmax(2L, as.integer(tbl))), names(tbl))
          if (length(vec) == 0L) stop("KmeansSolution2 returned no domain results.")
          kmeans_nstrata_vec(vec)
          kmeans_solution_stored(km)
          kmeans_sig(km_sig(method, maxcl))
        }, error = function(e) {
          showNotification(paste("KMeans2 error:", e$message), type = "error")
        })
      })
    } else {
      req(rv$strata, rv$errors)
      withProgress(message = "Running KmeansSolution for nStrata…", {
        tryCatch({
          km  <- KmeansSolution(strata = rv$strata, errors = rv$errors,
                                maxclusters = maxcl)
          tbl <- tapply(km$suggestions, km$domainvalue,
                        function(x) length(unique(x)))
          vec <- setNames(as.integer(pmax(2L, as.integer(tbl))), names(tbl))
          if (length(vec) == 0L) stop("KmeansSolution returned no domain results.")
          kmeans_nstrata_vec(vec)
          kmeans_solution_stored(km)
          kmeans_sig(km_sig(method, maxcl))
        }, error = function(e) {
          showNotification(paste("KMeans error:", e$message), type = "error")
        })
      })
    }
  }
  observeEvent(list(input$kmeans_nstrata, input$optim_method), run_kmeans())

  output$ui_kmeans_nstrata_info <- renderUI({
    if (!isTRUE(input$kmeans_nstrata)) return(NULL)
    vec <- kmeans_nstrata_vec()
    if (is.null(vec) || length(vec) == 0L) {
      helpText("Build strata (tab 3) first, then activate this option.")
    } else {
      helpText(paste("KMeans per domain:",
                     paste(seq_along(vec), vec, sep = ": ", collapse = ",  ")))
    }
  })

  observeEvent(input$btn_optim, {
    req(rv$strata, rv$errors, rv$frame)

    # Safety guard: continuous method is incompatible with categorised X variables
    if (!is.null(xfreq_df()) && input$optim_method == "continuous") {
      showNotification(
        "One or more X variables are categorised: only the 'atomic' method can be used. Please switch method.",
        type = "error", duration = 10)
      return()
    }

    optim_log("Running optimisation…\n")

    # Recompute the KMeans nStrata/seed if strata, errors, method or
    # maxclusters changed since it was last computed.
    if (isTRUE(input$kmeans_nstrata)) {
      cur_method <- if (is.null(input$optim_method)) "atomic" else input$optim_method
      cur_maxcl  <- max(2L, as.integer(if (is.null(input$nstrata)) 10L else input$nstrata))
      if (!identical(kmeans_sig(), km_sig(cur_method, cur_maxcl))) {
        optim_log(paste0(optim_log(),
          "KMeans solution outdated (strata, constraints or settings changed): recomputing...\n"))
        kmeans_nstrata_vec(NULL)
        kmeans_solution_stored(NULL)
        run_kmeans()
      }
    }

    withProgress(message = "Optimising stratification…", value = 0, {
      tryCatch({
        # Snapshot reactive values before the long computation.
        # When take-all strata are active: framesamp = frame minus framecens units
        # (must NOT contain any unit present in framecens), plus framecens is added.
        # Without take-all: framesamp = full frame, no framecens.
        use_cens_now <- isTRUE(input$use_cens) && !is.null(rv$framecens)
        frame_snap   <- if (use_cens_now) rv$framesamp else rv$frame
        cens_snap    <- if (use_cens_now) rv$framecens  else NULL
        errors_snap  <- rv$errors

        # errors has exactly one row per domain — more reliable than frame$domainvalue
        ndoms       <- nrow(errors_snap)

        # Determine nStrata vector and optional KMeans seed
        suggestions  <- NULL
        nstrata_user <- max(1L, as.integer(if (is.null(input$nstrata)) 5L else input$nstrata))

        use_km_vec <- isTRUE(input$kmeans_nstrata) && !is.null(kmeans_nstrata_vec())

        if (use_km_vec) {
          # KMeans ran on rv$strata — its result is the authoritative domain count.
          nstrata_vec <- kmeans_nstrata_vec()   # length == number of domains
          optim_log(paste0(optim_log(),
            "nStrata from KmeansSolution per domain: ",
            paste(names(nstrata_vec), nstrata_vec, sep = ": ", collapse = ", "),
            "\n"))

          if (isTRUE(input$use_kmeans)) {
            suggestions <- kmeans_solution_stored()
          }
        } else {
          dom_vals    <- as.character(errors_snap$domainvalue)
          nstrata_vec <- setNames(rep(nstrata_user, ndoms), dom_vals)
          optim_log(paste0(optim_log(), "nStrata (uniform): ", nstrata_user, "\n"))
          # Seeds (suggestions) are only supported by the atomic method.
          # For continuous, skip KmeansSolution entirely.
          if (isTRUE(input$use_kmeans) && input$optim_method == "atomic") {
            incProgress(0.1, detail = "KmeansSolution…")
            km_seed <- if (!is.null(kmeans_solution_stored())) {
              kmeans_solution_stored()
            } else {
              tryCatch(
                KmeansSolution(strata = rv$strata, errors = errors_snap,
                               maxclusters = nstrata_user),
                error = function(e) {
                  optim_log(paste0(optim_log(),
                    "KmeansSolution failed (continuing without seed): ", e$message, "\n"))
                  NULL
                }
              )
            }
            suggestions <- km_seed
          }
        }

        incProgress(0.2, detail = "Calling optimStrata…")

        if (input$optim_method == "atomic") {
          # Capture convergence plots via a multi-page PNG device.
          tmp_base <- file.path(tempdir(),
                                paste0("optim_plot_%03d_", as.integer(Sys.time()), ".png"))
          grDevices::png(tmp_base, width = 700, height = 380)
          sol <- tryCatch({
            # When take-all strata are active, disable the KMeans seed for the
            # atomic method: KmeansSolution is computed on the full strata/frame,
            # while atomic optimisation runs on framesamp (+ optional framecens).
            # Passing those suggestions can create length mismatches internally.
            if (use_cens_now && !is.null(suggestions)) {
              suggestions <- NULL
              optim_log(paste0(
                optim_log(),
                "Take-all strata active: disabling KMeans seed for atomic method.\n"
              ))
            }

            args_atomic <- list(
              method      = "atomic",
              errors      = errors_snap,
              framesamp   = frame_snap,
              iter        = input$iter,
              pops        = input$pops,
              nStrata     = nstrata_vec,
              suggestions = suggestions,
              showPlot    = TRUE,
              parallel    = FALSE   # strataGenalg not found in workers on 1.5-5
            )
            # Only add framecens when take-all units are actually defined —
            # passing framecens = NULL explicitly can break internal nrow() calls.
            if (!is.null(cens_snap)) args_atomic$framecens <- cens_snap
            do.call(optimStrata, args_atomic)
          }, error = function(e) {
            tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            stop(e)
          })
          tryCatch(grDevices::dev.off(), error = function(e) NULL)

          # rbga plots once per domain per iteration → keep only the last plot
          # per domain (every ndoms-th file from the end).
          generated <- sort_plot_files(Sys.glob(sub("%03d", "*", tmp_base, fixed = TRUE)))
          n_gen  <- length(generated)
          if (n_gen > 0 && n_gen >= ndoms && (n_gen %% ndoms) == 0) {
            step      <- n_gen %/% ndoms
            generated <- generated[seq(step, n_gen, by = step)]
          }
          rv$optim_plot_files <- if (length(generated) > 0) generated else NULL

        } else {
          # Continuous method — nStrata is a per-domain vector (not a scalar).
          # Any seed computed above (KmeansSolution) is not valid here: the only
          # supported seed is prepareSuggestion(), if the flag is ticked.
          suggestions <- NULL
          if (isTRUE(input$use_suggestions)) {
            incProgress(0.1, detail = "prepareSuggestion…")
            suggestions <- tryCatch(
              prepareSuggestion(frame = frame_snap, nStrata = nstrata_vec),
              error = function(e) {
                optim_log(paste0(optim_log(),
                  "prepareSuggestion failed (continuing without suggestions): ",
                  e$message, "\n"))
                NULL
              }
            )
            if (!is.null(suggestions))
              optim_log(paste0(optim_log(), "Using suggestions from prepareSuggestion().\n"))
          }
          tmp_base <- file.path(tempdir(),
                                paste0("optim_plot_%03d_", as.integer(Sys.time()), ".png"))
          grDevices::png(tmp_base, width = 700, height = 380)
          sol <- tryCatch({
            args_cont <- list(
              method    = "continuous",
              errors    = errors_snap,
              framesamp = frame_snap,
              iter      = input$iter,
              pops      = input$pops,
              nStrata   = nstrata_vec,
              showPlot  = TRUE,
              parallel  = FALSE
            )
            if (!is.null(cens_snap)) args_cont$framecens <- cens_snap
            if (!is.null(suggestions)) args_cont$suggestions <- suggestions
            optim_log(paste0(optim_log(),
              "framesamp rows: ", nrow(frame_snap),
              if (!is.null(cens_snap)) paste0("  |  framecens rows: ", nrow(cens_snap)) else "",
              "  |  nStrata: ", paste(names(nstrata_vec), nstrata_vec, sep="=", collapse=", "),
              "\n"))
            do.call(optimStrata, args_cont)
          }, error = function(e) {
            tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            stop(e)
          })
          tryCatch(grDevices::dev.off(), error = function(e) NULL)

          generated <- sort_plot_files(Sys.glob(sub("%03d", "*", tmp_base, fixed = TRUE)))
          n_gen  <- length(generated)
          if (n_gen > 0 && n_gen >= ndoms && (n_gen %% ndoms) == 0) {
            step      <- n_gen %/% ndoms
            generated <- generated[seq(step, n_gen, by = step)]
          }
          # The continuous method plots the domains in order of first appearance
          # in the frame (not sorted): reorder so that page k = domain k.
          if (length(generated) == ndoms) {
            dom_order <- unique(as.integer(as.character(
              c(frame_snap$domainvalue,
                if (!is.null(cens_snap)) cens_snap$domainvalue))))
            if (length(dom_order) == ndoms) generated <- generated[order(dom_order)]
          }
          rv$optim_plot_files <- if (length(generated) > 0) generated else NULL
        }
        rv$optim_plot_idx <- 1L

        rv$solution  <- sol

        # aggr_strata field name may differ by method/version; try both
        aggr <- if (!is.null(sol$aggr_strata)) sol$aggr_strata else sol$outstrata
        # Round SOLUZ to integers: selectSample passes SOLUZ to sample(size=...)
        # which truncates non-integers, so storing rounded integers ensures the
        # drawn sample size matches the displayed total.
        aggr$SOLUZ   <- round(aggr$SOLUZ)
        rv$outstrata <- aggr

        rv$framenew  <- sol$framenew   # frame with LABEL column needed by summaryStrata
        n_total <- sum(aggr$SOLUZ)

        used_suggestions <- !is.null(suggestions)
        opt_lines <- c(
          '# 7. Optimisation',
          paste0('nstrata <- ', named_numvec_to_r(nstrata_vec))
        )
        if (used_suggestions && identical(input$optim_method, "continuous")) {
          opt_lines <- c(opt_lines,
            paste0('suggestions <- prepareSuggestion(frame = ',
                   if (use_cens_now) 'framesamp' else 'frame',
                   ', nStrata = nstrata)')
          )
        } else if (used_suggestions) {
          if (isTRUE(input$kmeans_nstrata) && !is.null(kmeans_solution_stored())) {
            opt_lines <- c(opt_lines,
              paste0('kmeans_seed <- KmeansSolution(strata, errors, maxclusters = ',
                     nstrata_user, ')'),
              'suggestions <- kmeans_seed'
            )
          } else {
            opt_lines <- c(opt_lines,
              paste0('suggestions <- KmeansSolution(strata, errors, maxclusters = ',
                     nstrata_user, ')')
            )
          }
        }
        if (identical(input$optim_method, 'atomic')) {
          args_lines <- c(
            'args_optim <- list(',
            '  method = "atomic",',
            '  errors = errors,',
            if (use_cens_now) '  framesamp = framesamp,' else '  framesamp = frame,',
            paste0('  iter = ', input$iter, ','),
            paste0('  pops = ', input$pops, ','),
            '  nStrata = nstrata,',
            if (used_suggestions) '  suggestions = suggestions,' else NULL,
            '  showPlot = TRUE,',
            '  parallel = FALSE',
            if (use_cens_now) ',  framecens = framecens' else NULL,
            ')',
            'solution <- do.call(optimStrata, args_optim)'
          )
        } else {
          args_lines <- c(
            'args_optim <- list(',
            '  method = "continuous",',
            '  errors = errors,',
            if (use_cens_now) '  framesamp = framesamp,' else '  framesamp = frame,',
            paste0('  iter = ', input$iter, ','),
            paste0('  pops = ', input$pops, ','),
            '  nStrata = nstrata,',
            if (used_suggestions) '  suggestions = suggestions,' else NULL,
            '  showPlot = TRUE,',
            '  parallel = FALSE',
            if (use_cens_now) ',  framecens = framecens' else NULL,
            ')',
            'solution <- do.call(optimStrata, args_optim)'
          )
        }
        opt_lines <- c(opt_lines, args_lines,
          'outstrata <- if (!is.null(solution$aggr_strata)) solution$aggr_strata else solution$outstrata',
          'outstrata$SOLUZ <- round(outstrata$SOLUZ)',
          'framenew <- solution$framenew'
        )
        append_script_step(paste(opt_lines, collapse = '
'), '

')
        optim_log(paste0("Optimisation complete.\nTotal sample size: ", n_total,
                         "\n(Detailed domain output printed to R console)"))
        incProgress(1, detail = "Done.")
        showNotification("Optimisation complete!", type = "message")
      }, error = function(e) {
        # withCallingHandlers captures the stack before unwinding; tryCatch does not.
        # We re-run inside withCallingHandlers to get a real traceback.
        tb_env <- new.env(parent = emptyenv()); tb_env$calls <- ""
        withCallingHandlers(
          tryCatch(stop(e), error = function(e2) NULL),
          error = function(e2) {
            tb_env$calls <- paste(
              vapply(sys.calls(), function(x) deparse(x)[[1]], ""),
              collapse = "\n  ")
          }
        )
        optim_log(paste0("ERROR: ", e$message,
                         "\n\nCall stack:\n  ", tb_env$calls))
        showNotification(paste("Optimisation error:", e$message),
                         type = "error", duration = 15)
      })
    })
  })

  output$txt_optim_log <- renderPrint({ cat(optim_log()) })

  observeEvent(input$btn_plot_prev, {
    req(rv$optim_plot_files)
    rv$optim_plot_idx <- max(1L, rv$optim_plot_idx - 1L)
  })
  observeEvent(input$btn_plot_next, {
    req(rv$optim_plot_files)
    rv$optim_plot_idx <- min(length(rv$optim_plot_files), rv$optim_plot_idx + 1L)
  })

  output$txt_plot_page <- renderUI({
    req(rv$optim_plot_files)
    n <- length(rv$optim_plot_files)
    strong(paste("Domain", rv$optim_plot_idx, "of", n))
  })

  output$plot_optim <- renderImage({
    req(rv$optim_plot_files)
    idx <- rv$optim_plot_idx
    files <- rv$optim_plot_files
    idx <- max(1L, min(length(files), idx))
    list(src = files[[idx]], contentType = "image/png",
         width = "100%", alt = paste("Convergence plot domain", idx))
  }, deleteFile = FALSE)

  output$txt_expected_cv <- renderPrint({
    req(rv$outstrata)
    tryCatch(expected_CV(rv$outstrata), error = function(e) cat("Error:", e$message))
  })

  output$tbl_strata_summary <- renderDT({
    req(rv$framenew, rv$outstrata)
    tryCatch({
      ss <- summaryStrata(rv$framenew, rv$outstrata, progress = FALSE)
      datatable(ss, options = list(scrollX = TRUE, pageLength = 10))
    }, error = function(e) {
      datatable(data.frame(Error = e$message))
    })
  })

  # ── STEP 7: Results / Evaluation ─────────────────────────────────────────────
  output$ui_plot_vars <- renderUI({
    req(rv$framenew)
    # framenew uses uppercase column names (DOMAINVALUE, X1, X2, ...)
    xcols <- names(rv$framenew)[grep("^X[0-9]", names(rv$framenew))]
    fluidRow(
      column(6, selectInput("plot_x1", "X axis variable", xcols, xcols[1])),
      column(6, selectInput("plot_x2", "Y axis variable", xcols,
                            if (length(xcols) > 1) xcols[2] else xcols[1])),
      column(12, uiOutput("ui_plot_domain"))
    )
  })

  output$ui_plot_domain <- renderUI({
    req(rv$framenew)
    # framenew has uppercase DOMAINVALUE column
    dom_col <- if ("DOMAINVALUE" %in% names(rv$framenew)) "DOMAINVALUE" else "domainvalue"
    doms <- sort(unique(rv$framenew[[dom_col]]))
    selectInput("plot_domain", "Domain", choices = doms, selected = doms[1])
  })

  output$plot_strata2d <- renderPlot({
    req(rv$solution, rv$framenew, input$plot_x1, input$plot_x2)
    tryCatch({
      plotStrata2d(
        rv$framenew,
        rv$outstrata,
        domain = as.numeric(input$plot_domain),
        vars   = c(input$plot_x1, input$plot_x2)
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Plot error:", e$message), cex = 1.1)
    })
  })

  observeEvent(input$btn_eval, {
    req(rv$framenew, rv$outstrata)
    withProgress(message = "Evaluating solution…", {
      tryCatch({
        # ev <- evalSolution(rv$framenew, rv$outstrata,
        #                    nsampl = input$nsampl, progress = FALSE)
        ev <- evalSolution(
          rv$framenew,
          rv$outstrata,
          nsampl = input$nsampl,
          writeFiles = FALSE,
          progress = FALSE
        )
        rv$eval_res <- ev
        append_script_step('# Optional: evaluate solution
eval_res <- evalSolution(framenew, outstrata, nsampl = ', input$nsampl, ', progress = FALSE)
')
        showNotification("Evaluation complete.", type = "message")
      }, error = function(e) {
        showNotification(paste("Eval error:", e$message), type = "error", duration = 10)
      })
    })
  })

  output$plot_eval_cv <- renderPlot({
    req(rv$eval_res)
    # coeff_var: columns CV1, CV2, ..., dom (one row per domain)
    cv      <- rv$eval_res$coeff_var
    cv_cols <- setdiff(names(cv), "dom")

    df <- do.call(rbind, lapply(cv_cols, function(v) {
      data.frame(domain   = as.character(cv$dom),
                 variable = v,
                 value    = cv[[v]],
                 stringsAsFactors = FALSE)
    }))

    # One threshold per domain and CV variable, taken from rv$errors
    # (columns CV1, CV2, ... and domainvalue; constraints may differ by domain).
    thresh <- NULL
    if (!is.null(rv$errors) && "domainvalue" %in% names(rv$errors)) {
      thresh <- do.call(rbind, lapply(intersect(cv_cols, names(rv$errors)), function(v) {
        # evalSolution labels domains "DOM1", "DOM2", ...: use the same labels
        dom_raw   <- as.character(rv$errors$domainvalue)
        dom_label <- if (all(paste0("DOM", dom_raw) %in% as.character(cv$dom)))
                       paste0("DOM", dom_raw) else dom_raw
        data.frame(domain    = dom_label,
                   variable  = v,
                   threshold = rv$errors[[v]],
                   stringsAsFactors = FALSE)
      }))
      if (!is.null(thresh))
        thresh <- thresh[!is.na(thresh$threshold), , drop = FALSE]
    }

    p <- ggplot(df, aes(x = domain, y = value, fill = variable)) +
      geom_col(position = "dodge") +
      labs(title = "CV by domain and variable  (black ticks = constraints)",
           x = "Domain", y = "CV", fill = "Variable") +
      theme_bw()

    if (!is.null(thresh) && nrow(thresh) > 0) {
      # black tick over each bar, dodged exactly like the bars
      p <- p + geom_errorbar(data = thresh,
                             aes(x = domain, ymin = threshold, ymax = threshold,
                                 group = variable),
                             position = position_dodge(width = 0.9),
                             width = 0.9, linewidth = 0.8, colour = "black",
                             inherit.aes = FALSE)
    }
    p
  })

  output$plot_eval_bias <- renderPlot({
    req(rv$eval_res)
    # rel_bias: columns y1, y2, ..., dom (one row per domain)
    bias     <- rv$eval_res$rel_bias
    bias_cols <- setdiff(names(bias), "dom")
    df <- do.call(rbind, lapply(bias_cols, function(v) {
      data.frame(domain   = as.character(bias$dom),
                 variable = v,
                 value    = bias[[v]],
                 stringsAsFactors = FALSE)
    }))
    ggplot(df, aes(x = domain, y = value, fill = variable)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Relative Bias by domain and variable", x = "Domain", y = "Rel. Bias",
           fill = "Variable") +
      theme_bw()
  })

  observeEvent(input$btn_adjust, {
    req(rv$outstrata)
    target <- input$adj_size
    if (target <= 0) {
      showNotification("Enter a positive target sample size.", type = "warning")
      return()
    }
    tryCatch({
      adj_out  <- capture.output(
        new_strata <- adjustSize(size = target, strata = rv$outstrata)
      )
      # Round SOLUZ to integers; then correct rounding accumulation by
      # adjusting ±1 per stratum to reach the target. Take-all strata
      # (CENS == 1) are never modified: adjustSize() already accounts for them.
      new_strata$SOLUZ <- round(new_strata$SOLUZ)
      is_cens  <- if ("CENS" %in% names(new_strata)) new_strata$CENS == 1 else
                    rep(FALSE, nrow(new_strata))
      actual_size <- sum(new_strata$SOLUZ)
      delta <- actual_size - target
      if (delta > 0) {
        # Too many: reduce strata with largest SOLUZ by 1, keeping minimum 2
        eligible <- which(new_strata$SOLUZ > 2 & !is_cens)
        if (length(eligible) > 0) {
          reduce_n   <- min(delta, length(eligible))
          reduce_idx <- eligible[order(-new_strata$SOLUZ[eligible])][seq_len(reduce_n)]
          new_strata$SOLUZ[reduce_idx] <- new_strata$SOLUZ[reduce_idx] - 1
        }
      } else if (delta < 0) {
        # Too few: increase strata with smallest SOLUZ by 1 (only where N allows)
        eligible <- which(!is_cens)
        if ("N" %in% names(new_strata))
          eligible <- eligible[new_strata$SOLUZ[eligible] < new_strata$N[eligible]]
        if (length(eligible) > 0) {
          increase_n   <- min(-delta, length(eligible))
          increase_idx <- eligible[order(new_strata$SOLUZ[eligible])][seq_len(increase_n)]
          new_strata$SOLUZ[increase_idx] <- new_strata$SOLUZ[increase_idx] + 1
        }
      }
      # Cap SOLUZ at stratum population N (can't sample more than the population)
      if ("N" %in% names(new_strata))
        new_strata$SOLUZ <- pmin(new_strata$SOLUZ, new_strata$N)
      actual_size <- sum(new_strata$SOLUZ)
      rv$outstrata <- new_strata
      append_script_step('# Optional: adjust sample size
outstrata <- adjustSize(size = ', target, ', strata = outstrata)
outstrata$SOLUZ <- round(outstrata$SOLUZ)
is_cens <- if ("CENS" %in% names(outstrata)) outstrata$CENS == 1 else rep(FALSE, nrow(outstrata))
delta <- sum(outstrata$SOLUZ) - ', target, '
if (delta > 0) {
  eligible <- which(outstrata$SOLUZ > 2 & !is_cens)
  if (length(eligible) > 0) {
    idx <- eligible[order(-outstrata$SOLUZ[eligible])][seq_len(min(delta, length(eligible)))]
    outstrata$SOLUZ[idx] <- outstrata$SOLUZ[idx] - 1
  }
} else if (delta < 0) {
  eligible <- which(!is_cens)
  if ("N" %in% names(outstrata)) eligible <- eligible[outstrata$SOLUZ[eligible] < outstrata$N[eligible]]
  if (length(eligible) > 0) {
    idx <- eligible[order(outstrata$SOLUZ[eligible])][seq_len(min(-delta, length(eligible)))]
    outstrata$SOLUZ[idx] <- outstrata$SOLUZ[idx] + 1
  }
}
if ("N" %in% names(outstrata)) outstrata$SOLUZ <- pmin(outstrata$SOLUZ, outstrata$N)
')
      # Invalidate previous evaluation so re-running Simulation uses updated SOLUZ
      rv$eval_res  <- NULL
      msg <- paste0("Size adjusted. Requested: ", target,
                    " | Achieved: ", actual_size,
                    if (actual_size != target)
                      " (minimum 2 per stratum or stratum population size limits total)"
                    else "")
      showNotification(msg, type = "message", duration = 8)
    }, error = function(e) {
      showNotification(paste("Adjust error:", e$message), type = "error", duration = 10)
    })
  })

  output$tbl_adj_strata <- renderDT({
    req(rv$outstrata)
    datatable(rv$outstrata, options = list(scrollX = TRUE, pageLength = 10))
  })

  # ── STEP 8: Sample selection ─────────────────────────────────────────────────
  output$ui_sort_var <- renderUI({
    req(rv$framenew)
    cols <- names(rv$framenew)
    conditionalPanel("input.sample_method == 'systematic'",
      selectInput("sort_var", "Sort variable (systematic)", cols, cols[1])
    )
  })

  observeEvent(input$btn_select, {
    req(rv$framenew, rv$outstrata)
    tryCatch({
      if (input$sample_method == "random") {
        samp <- selectSample(rv$framenew, rv$outstrata, writeFiles = FALSE)
      } else {
        samp <- selectSampleSystematic(rv$framenew, rv$outstrata,
                                       sortvariable = input$sort_var)
      }
      rv$sample_df <- samp
      if (input$sample_method == "random") {
        append_script_step('# 8. Select sample
sample_df <- selectSample(framenew, outstrata, writeFiles = FALSE)
')
      } else {
        append_script_step('# 8. Select sample
sample_df <- selectSampleSystematic(framenew, outstrata, sortvariable = ', quote_r_string(input$sort_var), ')
')
      }
      showNotification(paste("Sample selected: n =", nrow(samp)), type = "message")
    }, error = function(e) {
      showNotification(paste("Sample error:", e$message), type = "error", duration = 10)
    })
  })

  output$tbl_sample <- renderDT({
    req(rv$sample_df)
    datatable(head(rv$sample_df, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$txt_sample_summary <- renderPrint({
    req(rv$sample_df)
    cat("Total sample size:", nrow(rv$sample_df), "\n\n")
    print(table(rv$sample_df$DOMAINVALUE))
  })

  build_script_lines <- reactive({
    steps <- rv$script_steps
    if (is.null(steps) || length(steps) == 0) {
      return(c(
        '# No executable workflow has been recorded yet.',
        '# Use the app buttons first, then download the script again.'
      ))
    }
    c(
      '# SamplingStrata workflow script generated by the Shiny app.',
      '#',
      '# What this script does:',
      '# - loads the main input dataset',
      '# - prepares the selected variables',
      '# - builds the SamplingStrata frame and strata objects',
      '# - defines precision constraints',
      '# - optionally loads and validates a take-all subset, then builds framecens',
      '# - runs the optimisation and stores framenew and outstrata',
      '# - optionally evaluates the solution and/or adjusts the sample size',
      '# - selects the final sample',
      '#',
      '# Expected input files:',
      '# - the main dataset specified in the data-loading section below',
      '# - optionally, a take-all subset CSV file if the take-all step was used',
      '#',
      '# Output objects/files:',
      '# - framenew: frame with stratum labels returned by optimStrata()',
      '# - outstrata: optimised strata allocation table',
      '# - sample_df: final selected sample',
      '# - optional file export at the end: sample.csv',
      '',
      steps,
      '',
      '# Optional export of the selected sample',
      'if (exists("sample_df")) write.csv(sample_df, "sample.csv", row.names = FALSE)'
    )
  })

  output$txt_script_preview <- renderText({
    paste(build_script_lines(), collapse = "
")
  })

  output$btn_download <- downloadHandler(
    filename = "sample.csv",
    content  = function(file) {
      req(rv$sample_df)
      write.csv(rv$sample_df, file, row.names = FALSE)
    }
  )

  output$btn_download_script <- downloadHandler(
    filename = function() {
      paste0("samplingstrata_workflow_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    },
    content = function(file) {
      writeLines(build_script_lines(), con = file, useBytes = TRUE)
    }
  )

}

shinyApp(ui = ui, server = server)
