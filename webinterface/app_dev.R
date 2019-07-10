#!/usr/bin/Rscript
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(shinyWidgets)
library(dashboardthemes) # library(devtools); install_github("nik01010/dashboardthemes")
library(shinyhelper)
library(plotly)
library(viridis)

# welcome
welcome_page <- dashboardBody(
  fluidRow(
 #   box(title = h2("Overview"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
    boxPlus(title = h2("Overview"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
        h3(strong("fastABC"), "is a DNA sequence analysis workflow to study the demographic history of sampled populations or species by using Approximate Bayesian Computations."),
        h3("From a single uploaded input file containing sequenced genes or DNA fragments,", strong("fastABC"), "will:"),
        h3(strong("1."), "simulate different models/scenarios."),
        h3(strong("2."), "select the best model using an ABC approach based on", a(span(strong("random forests."), style = "color:teal"), href="https://cran.r-project.org/web/packages/abcrf/index.html", target="_blank")),
        h3(strong("3."), "estimate the parameters of the best model using a", a(span(strong("neural network"), style = "color:teal"), href="https://cran.r-project.org/web/packages/abc/index.html", target="_blank"), "approach."),
        h3(strong("4."), "measure the robustness of the analyses.", strong("fastABC"), "is transparent on the ability of its inferences to reproduce the observed data."),
        hr(),
        h3("The first goal of", strong("fastABC"), "is to distinguish between isolation versus migration models for sister gene pools."),
        h3("Its ultimate goal is to produce for each studied gene the probability of being associated with a species barrier.")
        
      )
    ),
  
  fluidRow(
    #box(title = h2("Compared demographic models"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
    boxPlus(title = h2("Compared demographic models"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                mainPanel(width=NULL, htmlOutput("models_picture"),
                         h3(strong("2 populations/species")),
                         h3(strong("SI"), "= strict isolation: subdivision of an ancestral diploid panmictic population (of size Nanc) in two diploid populations (of constant sizes Npop1 and Npop2) at time Tsplit."),
                         h3(strong("AM"), "= ancestral migration: the two newly formed populations continue to exchange alleles until time TAM."),
                         h3(strong("IM"), "= isolation with migration: the two daughter populations continuously exchange alleles until present time."),
                         h3(strong("SC"), "= secondary contact: the daughter populations first evolve in isolation (forward in time), then experience a secondary contact and start exchanging alleles at time TSC. Red phylogenies represent possible gene trees under each alternative model."),
                         hr(),
                         h3(strong("4 populations/species")),
                         h3("A single generalist model, declined in 64 sub-models according to if there is one:"),
                         h3("migration (bidirectional) or no migration between A and B, and/or between C and D."),
                         h3("bidirectional, unidirectional or no migration between A and C, and/or between B and D."),
                         h3("In case of migration between A and B (and between C and D), the gene flow takes place since their separation Tsplit_AB (and Tsplit_CD)."),
                         h3("In case of migration between A and C (and between B and D), the gene flow occurs during a secondary contact TSC_AC (and TSC_BD) lower than", code("min(c(Tsplit_AB, Tsplit_CD))"), ".")
                )
      
    ),
    
    
    #box(title = h2("Compared genomic models"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
    boxPlus(title = h2("Compared genomic models"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                
      fluidRow(
        # Sidebar panel for inputs
        column(width=4,
               # Input: Slider for the number of bins
               sliderInput(inputId = "Ne", label = h3("Effective population size:"), min = 0, max = 1000000, value = 10000, step=1000),
               sliderInput(inputId = "alpha", label = h3("Shape parameter #1:"), min = 1, max = 10, value = 10, step=0.1),
               sliderInput(inputId = "beta", label = h3("Shape parameter #2:"), min = 1, max = 10, value = 3, step=0.1)	
        ),
        
        # Main panel for displaying outputs
        column(width=8,
               # Output: density
               plotOutput(outputId = "genomic_hetero")
        )
      ),
      
      fluidRow(
        column(width=12,
               h3("All demographic models exist under AT LEAST two alternative genomic models:"),
               h3(strong("1."), "a model where the effective size", strong("(Ne)"), "is genomically homogeneous (orange bar), i.e., all locus are simulated by sharing the same", strong("Ne"), "value. In this model,", strong("fastABC"), "will try to estimate the value of", strong("Ne"), "best explaining the observed data.", strong("Ne"), "being independent between all populations (current, past)."),
               h3(strong("2."), "a model where", strong("Ne"), "is genomically heterogeneous (green distribution), i.e., all locus are simulated with a value of", strong("Ne"), "drawn in a Beta distribution. In this model,", strong("fastABC"), "will try to estimate the value of", strong("Ne"), "as well as the two shape parameters", strong("(shape1 and shape2)"), "that best explain the observations. Here,", strong("fastABC"), "assumes that all populations (current and past)", strong("share the same Beta distribution"), "but are independently rescaled by different", strong("Ne"), "values."),
               hr(),
               h3("In addition, all demographic models with migration have two alternative models of introgression:"),
               h3(strong("1."), "a model where all of the loci share the same introgression rate for a given direction, but these rates are independent between directions. Here,", strong("fastABC"), "will simply try to estimate the introgression rate for each direction."),
               h3(strong("2."), "a model where introgression rates are Beta distributed throughout genomes.", strong("fastABC"), "assumes independent Beta distributions for each direction where gene flow occurs."),
               htmlOutput("homo_hetero")
               
        )     
        )
    ),
    
#    box(title = h2("Model comparisons for 2 populations/species"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
    boxPlus(title = h2("Model comparisons for 2 populations/species"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
          htmlOutput("model_comparisons"),
          h3(strong("fastABC"), "performs hierarchical model comparisons."),
          h3(strong("1."), "comparison between all models with", strong("current isolation"), "({SI; AM} x {Ne_homo; Ne_hetero} x {M_homo; M_hetero}) versus", strong("ongoing migration"), "({IM; SC} x {Ne_homo; Ne_hetero} x {M_homo; M_hetero})"),
          h3(strong("2. if current isolation ->"), "comparison between", strong("SI"), "({Ne_homo; Ne_hetero}) versus", strong("AM"), "({Ne_homo; Ne_hetero} x {M_homo; M_hetero})"),
          h3(strong("2. if ongoing migration ->"), "comparison between", strong("IM"), "({Ne_homo; Ne_hetero} x {M_homo; M_hetero}) versus", strong("SC"), "({Ne_homo; Ne_hetero} x {M_homo; M_hetero})"),
          h3(strong("3."), "the last step is to determine whether effective size", strong("(Ne)"), "and migration rates", strong("(N.m)"), "are homogeneously or heterogenously distributed in genomes.")
          
    ),

    boxPlus(title = h2("Architecture"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
        column(width=12,
          h3("FastABC is composed of two elements:"),
          h3(strong("1."), "A web interface developed in", a(span(strong("Shiny,"), style = "color:teal"), href="https://shiny.rstudio.com/", target="_blank"), "which will execute..."),
          h3(strong("2."), "...a workflow managed by", a(span(strong("Snakemake."), style = "color:teal"), href="https://snakemake.readthedocs.io/en/stable/", target="_blank")),
          htmlOutput("welcome_picture"),
          hr(),
          h3("The code is fully open-source, freely distributed on", a(span(strong("GitHub"), style = "color:teal"), href="https://github.com/popgenomics/ABConline", target="_blank"), "and can be immediately redeployed on any cluster using", a(span(strong("SLURM"), style = "color:teal"), href="https://slurm.schedmd.com/documentation.html", target="_blank"), "thanks to a", a(span(strong("Singularity"), style = "color:teal"), href="https://sylabs.io/docs/#doc-3.2", target="_blank"), "image."),
          
          h3("This redeployment allows the user to modify the models to be compared, to add summary statistics, etc..."),
          h3("The workflow can be simply executed from the command line without going through the web interface."),
          
          h3("However, if desired, the web interface can also be freely hosted and linked to any cluster.")
        )
      )
  )
)

# upload
upload_data <- dashboardBody(
        tags$head(tags$style(HTML("a {color: black}"))),
        fluidRow(
         #   box(title = h2("Number of ABC analysis to run"), height = 250, width = 6, solidHeader = TRUE, background = NULL, status = "primary",
            boxPlus(title = h2("Number of ABC analysis to run"), height = 250, width = 6, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                shinyjs::useShinyjs(),
                #numericInput(inputId = "number_of_ABC", label = NULL,  value = 1, min = 1, max = 5, step = 1),
                selectInput("number_of_ABC", label = h4("1 to 5 ABC analyses can be performed from the same input file"), choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5), selected = 1)
                
               # actionButton("number_of_ABC_validation", "Validate")
            ),
            
        #    box(title = h2("Email address"), height = 250,  width = 6, solidHeader = TRUE, status = "primary",
            boxPlus(title = h2("Email address"), height = 250, width = 6, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                textInput("mail_address", label = h4("address to receive the download link of the results"), value = "user@gmail.com")
            )
        ),

     fluidRow(NULL, soldHeader = TRUE, status ="danger",
#              box(title = h2("Sequence Alignment Upload"), height = 200,  width = 6, solidHeader = TRUE, status = "success",
              boxPlus(title = h2("Sequence Alignment Upload"), height = 200,  width = 6, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,

                  fileInput("infile", label = NULL)
              ),
              
             
#             box(title = h2("Genomic regions"), height = 200,  width = 6, solidHeader = TRUE, status = "warning",
             boxPlus(title = h2("Genomic regions"), height = 200,  width = 6, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                 #                radioButtons("region", label = NULL, selected = "noncoding", choices = list("coding" = "coding", "non coding" = "noncoding"))
                 prettyRadioButtons("region", label = NULL, shape = "round", status = "warning", fill = TRUE, inline = TRUE, animation = "pulse", bigger = TRUE,
                                    selected = "noncoding", choices = list("coding" = "coding", "non coding" = "noncoding"))
                 )
        ),     
        
        fluidRow(align="left",
            #box(title = h2("Input file"), width = 6, solidHeader = TRUE, background = "green", status = "success",
                boxPlus(title = h2("Expected input file's format"), width = 6, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                h3(strong("Fasta file")),
                h3("A single fasta file containing all sequences obtained from all populations/species, and for all genes is the only inputfile to upload."),
                h3("Even sequences obtained from non-studied species can be included in the file. The user will specify the names of the species to consider after the upload ."),
                h3("Its format is largely inspired by the output of", a(strong("Reads2snp"), href="https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools", target="_blank"), " but can be post-produced without Reads2snp."),
                hr(),
                h3("All sequences's have to respect the following structure:"),
                h3(strong(">gene|species or population|individual|allele1 or allele2")),
                h3(strong("GTGATGCGTGTAGTCATG")),
                h3("With missing data only encoded by 'N'"),
                br()
            ),
            
            #box(title = h2("Informations about the uploaded file"),  width = 6, solidHeader = TRUE, status = "success",
            boxPlus(title = h2("Informations about the uploaded file"), width = 6, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                uiOutput("upload")
                
            )
        ),

        fluidRow(
            box("", width = 12, solidHeader = TRUE, status = "info",
                prettyCheckbox(inputId = "check_upload", shape = "round", value = FALSE,
                           label = strong("Please check/valid your choices"), icon = icon("check"),
                           animation = "tada", status = "success", bigger = TRUE)
            )
        ),

        fluidRow(
            boxPlus(
                title = h2("Example"), width = 12, icon = NULL, solidHeader = TRUE, background = NULL,
                boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                enable_label = TRUE, label_text = "CLICK TO DISPLAY AN EXAMPLE OF INPUT FILE", label_status = "success",
                p(">Hmel210004_196|chi|chi.CJ560|allele1"),
                p("NNNNNNNGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|chi|chi.CJ560|allele2"),
                p("NNNNNNNGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|chi|chi.CJ564|allele1"),
                p("NTGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|chi|chi.CJ564|allele2"),
                p("NTGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATATTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|flo|flo.CS2338|allele1"),
                p("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
                p(">Hmel210004_196|flo|flo.CS2338|allele2"),
                p("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
                p(">Hmel210004_196|flo|flo.CS2341|allele1"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|ros|ros.CJ2071|allele1"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|ros|ros.CJ2071|allele2"),
                p("ATGTCTCGGCCAGTGTTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCTGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|ros|ros.CJ531|allele1"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCTGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|num|nu_sil.MJ09-4125|allele1"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTGGAATTATACAAA"),
                p(">Hmel210004_196|num|nu_sil.MJ09-4125|allele2"),
                p("ATGTCTCGGCCAGTATTATCTACACACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTGGAATTATACAAA"),
                p(">Hmel210004_196|num|nu_sil.MJ09-4184|allele1"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
                p(">Hmel210004_196|num|nu_sil.MJ09-4184|allele2"),
                p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAATGGAATTTTCGTGGAATTATACAAA"),
                p(">Hmel219015_26|chi|chi.CAM25091|allele1"),
                p("GGAAATNNAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTNTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACANGAGAAAATAA"),
                p(">Hmel219015_26|chi|chi.CAM25091|allele2"),
                p("GGAAATNNAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTNTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTAAAATTTCACANGAGAAAATAA"),
                p(">Hmel219015_26|chi|chi.CAM25137|allele1"),
                p("GGAAATNNAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCANAGGAGAAAATAA"),
                p(">Hmel219015_26|chi|chi.CAM25137|allele2"),
                p("GGAAATNNAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTAAAATTTCANAAGAGAAAATAA"),
                p(">Hmel219015_26|flo|flo.CS12|allele1"),
                p("GGAAATGAAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAAGAGAAAATAA"),
                p(">Hmel219015_26|flo|flo.CS12|allele2"),
                p("GGAAATGAAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGCCTAGAAGCTGGAACAAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAAGAGAAAATAA"),
                p(">Hmel219015_26|flo|flo.CS13|allele1"),
                p("NNNNNNGAAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
                p(">Hmel219015_26|ros|ros.CAM1841|allele1"),
                p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
                p(">Hmel219015_26|ros|ros.CAM1841|allele2"),
                p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
                p(">Hmel219015_26|ros|ros.CAM1880|allele1"),
                p("NNNNNNNNNNNNNNNNNTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
                p(">Hmel219015_26|ros|ros.CAM1841|allele1"),
                p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
                p(">Hmel219015_26|ros|ros.CAM1841|allele2"),
                p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
                p(">Hmel219015_26|ros|ros.CAM1880|allele1"),
                p("NNNNNNNNNNNNNNNNNTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
                p(">Hmel219015_26|num|nu_sil.MJ09-4125|allele1"),
                p("NNNNNNGNNNNCTTTTGTATNAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTNAAATTTCACAAGAGAAAATAA"),
                p(">Hmel219015_26|num|nu_sil.MJ09-4125|allele2"),
                p("NNNNNNTNNNNCTTTTGTATNAATTCTGTTGAGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGACCTGGTCTTTCGCACTGATATTATATTACGAACTATTGGACAACCAGTGTACGTNAAATTTCACAAAAGAAAATAA"),
                p(">Hmel219015_26|num|nu_sil.MJ09-4184|allele1"),
                p("NGANATGAAAACNTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACNAAGCCATCTGATCTGGTCTTCCGCACTGATATTNTATTGCGAACTATGGGACAACCAATTTACGTNAAATTTCACANNAGAAAATAA"),
                p(">Hmel219015_26|num|nu_sil.MJ09-4184|allele2"),
                p("NGANATTAAAACNTTTGTATCAATTCTGTTGAGGCGATTTCGTCTAGAAGCTGTAACNAAGCCATCTGATCTGGTCTTTCGCACTGATATTNTATTACGAACTATTGGACAACCAGTGTACGTNAAATTTCACANNAGAAAATAA"),
                p("etc ...")
            ),
            
            boxPlus(
                title = h2("Description of the example"), width = 12, icon = NULL, solidHeader = TRUE, gradientColor = "success",
                boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                enable_label = TRUE, label_text = "CLICK TO DISPLAY THE DESCRIPTION OF EXPECTED INPUT FILE", label_status = "success",
                h3("Two genes are displayed from this file, they are named: ", strong("Hmel210004_196"), " and", strong("Hmel219015_26.")),
                br(),
                h3("Four populations are present in this file, named: ", strong("chi, flo, ros and num.")),
                h3("Only species whose names are specified in the ", strong("Populations/species"), " menu are considered. This does not prevent the uploaded file from containing other species."),
                br(),
                h3("Two diploid individuals are sequenced by species/population. This number is obviously allowed to vary between species/populations, according to the sequencing strategy and its success.")
                
            )
        )
)


filtering <- dashboardBody(
    fluidRow(
        column(width = 4,
               #box(title = h2("Maximum proportion of N"), width = NULL, solidHeader = TRUE, status = "primary", height = 225,
                boxPlus(title = h2("Maximum proportion of missing data (N, gaps, ...)"), height = 225, width = NULL, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                   sliderInput("max_N_tolerated", label = NULL,  min = 0, max = 1, value = 0.1, step = 0.005)
               ),
               boxPlus(
                   title = h3("max_N_tolerated"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "EXPLANATIONS", label_status = "primary",
                    h3("-Variable between 0 and 1."),
                    h3("-Defines the maximum proportion of N in the sequence of a gene in an individual beyond which this sequence is not considered."),
                    hr(),
                    h3(a(span(strong("Example", style = "color:blue")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/max_N_tolerated.png", target="_blank"))
               )
        ),
        
        column(width = 4,
              # box(title = h2("Minimum gene length"), width = NULL, solidHeader = TRUE, status = "success", height = 225,
              boxPlus(title = h2("Minimum sequence length per gene"), height = 225, width = NULL, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                   numericInput("Lmin", label = NULL, value = 30)
               ),

               boxPlus(
                    title = h3("Lmin"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
                    boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                    enable_label = TRUE, label_text = "EXPLANATIONS", label_status = "success",
                   
                    h3("-Positive integer (>0)"),
                    h3("-Minimum number of treatable sites below which a gene is removed from the analysis."),
                    br(),
                    h3("In a noncoding sequence: a site is an alignment of nucleotides for a single given nucleotide position, including all individuals among the species considered."),
                    br(),
                    h3("In a coding sequence: a site is an alignment for a given codon, comprising all the individuals among the species considered."),
                    h3("A coding position is not considered if:"),
                    h3("    -a codon alignment contains a non-synonymous polymorphism."),
                    h3("    -more than two codons segregate (even synonyms)."),
                    h3("    -at least one N is found in a codon, in an individual."),
                    br(),
                    h3("Number of positions to consider = (number of ", strong("monomorphic positions"), "that can be considered) + (number of ", strong("biallelic positions"), "that can be considered)."),
                    hr(),
                    h3(a(span(strong("Example", style = "color:green")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/Lmin.png", target="_blank"))
               )
        ),
        
        column(width = 4,
               #box(title = h2("Minimum number of sequences"), width = NULL, solidHeader = TRUE, status = "warning", height = 225,
                boxPlus(title = h2("Minimum number of sequences per gene and per population/species"), height = 225, width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                   numericInput("nMin", label = NULL, value = 12)
               ),
               boxPlus(
                   title = h3("nMin"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "EXPLANATIONS", label_status = "warning",
                   h3(strong("fastABC"), " starts for each gene by eliminating individual sequences containing too many N and gaps", span(strong("(max_N_tolerated; blue box)", style = "color:blue")), "."),
                   br(),
                   h3("If for a gene and", strong("within a population/species"), "there are fewer ", strong("nMin"), "sequences left, then the gene is not considered in the ", strong("ABC"), "analysis."),
                   br(),
                   h3(strong("If an outgroup is specified:")),
                   h3(strong("nMin"), " becomes the number of sequences sampled for each species, for each locus, to produce a standardized joint SFS used by ", strong("ABC"), "."),
                   hr(),
                   h3(a(span(strong("Example", style = "color:orange")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/nMin.png", target="_blank"))
               )
        )
    ),
    
    fluidRow(
        box("", width = 12, solidHeader = TRUE, status = "info",
            prettyCheckbox(inputId = "check_filtering", shape = "round", value = FALSE,
                           label = strong("Please check/valid your choices"), icon = icon("check"),
                           animation = "tada", status = "success", bigger = TRUE)
        )
    )
)


populations <- dashboardBody(
    fluidRow(
        #box(title = h2("Number of populations/species"), width = 4, solidHeader = TRUE, status = "primary", height = 600,
        boxPlus(title = h2("Number of populations/species"), height = NULL, width = 6, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                    
#            radioButtons("nspecies", label = h3("Number of gene pools"), choices = list("One gene pool" = 1, "Two gene pools" = 2, "Four gene pools" = 4), selected = 2),
            prettyRadioButtons("nspecies", label = h3("Number of gene pools"), shape = "round", status = "primary", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
                               choices = list("One gene pool" = 1, "Two gene pools" = 2, "Four gene pools" = 4), selected = 2),
            uiOutput("input_names_ui")
        ),
        
        #box(title = h2("Outgroup species"), width = 4, solidHeader = TRUE, status = "warning", height = 600,
        boxPlus(title = h2("Outgroup species"), height = NULL, width = 6, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                    
#            radioButtons("presence_outgroup", label = h3("Presence of an outgroup"), choices = list("no" = "no", "yes" = "yes"), selected = "no"),
            prettyRadioButtons("presence_outgroup", label = h3("Presence of an outgroup"), shape = "round", status = "warning", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
                               choices = list("no" = "no", "yes" = "yes"), selected = "no"),
            uiOutput("input_names_outgroup_ui")
        )
    ),

    fluidRow(
        box("", width = 12, solidHeader = TRUE, status = "info",
            prettyCheckbox(inputId = "check_populations", shape = "round", value = FALSE,
                           label = strong("Please check/valid your choices"), icon = icon("check"),
                           animation = "tada", status = "success", bigger = TRUE)
        )
    )
)


prior <- dashboardBody(
    fluidRow(
        column(width = 6,
              #box(title = h2("Mutation and recombination"), width = NULL, solidHeader = TRUE, status = "primary", height = 250,
              boxPlus(title = h2("Mutation and recombination"), height = NULL, width = NULL, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                   fluidRow(
                        column(width=5, numericInput("mu", label = h5('Mutation rate'), value = 0.000000003)),
                        column(width=5, numericInput("rho_over_theta", label = h5('Ratio r/µ'), value = 0.1))
                   )                ),
               
               boxPlus(
                   title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "teal",
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "Informations about mutation and recombination", label_status = "primary",
                   
                   h3("The mutation rate ", strong("(µ)"), "is", strong("the probability per generation and per nucleotide"), "that an allele will not be properly replicated."),
                   h3("If an external group ", strong("is not specified"), "then all genes/contigs/locus share the same µ."),
                   h3("If an external group ", strong("is specified"), "then the local µ, for a locus ", em("i"), " is corrected by ", strong("µ * div_i / div_avg"), " where ", strong("div_i"), " is the local divergence between the ingroup and the outgroup at that locus, and ", strong("div_avg"), " is the divergence averaged over loci"),
                   br(),
                   h3("The r/µ ratio is the ratio of recombination (/bp /generation) over mutation (/bp /generation)."),
                   h3("If the ratio is (unnecessarily) setted to values above 10, simulations will take too long.")
               ),
               
#               box(title = h2("Population sizes"), width = NULL, solidHeader = TRUE, status = "danger", height = 250,
               boxPlus(title = h2("Population size"), height = 250, width = NULL, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                    fluidRow(
                        column(width=5, numericInput("N_min", label = h5('min'), value = 100)),
                        column(width=5, numericInput("N_max", label = h5('max'), value = 1000000))
                   )
               ),
               
               boxPlus(
                   title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "danger",
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "Informations about population size", label_status = "danger",
                   
                   h3("The effective population size ", em(strong("Ne")), "is the number of diploid individuals within current and ancestral species/populations."),
                   hr(),
                   h3("In the", strong("ABC"), "simulations,", em(strong("Ne")), "will be drawn from the setted prior distribution independently for all current and ancestral species/populations")
               )
        ),
        
        column(width = 6,
               #box(title = h2("Time of split"), width = NULL, solidHeader = TRUE, status = "warning", height = 250,
              boxPlus(title = h2("Time of split"), height = NULL, width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                           
                  fluidRow(
                       column(width=5, numericInput("Tsplit_min", label = h5('min'), value = 100)),
                       column(width=5, numericInput("Tsplit_max", label = h5('max'), value = 1000000))
                   )
               ),

               boxPlus(
                   title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "warning",
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "Informations about the time of split", label_status = "warning",
                   
                   h3("The speciation time", strong(em("Tsplit")), "is expressed", strong("in number of generations.")),
                   h3("For annual organisms: one generation = one year."),
                   h3("For perennial organisms: one generation = average age for an individual to transmit a descendant (which is different from the age of sexual maturity)."),
                   br(),
                   h3("The ", em("prior"), " distribution is uniform between", strong(em("Tsplit_min")), "and", strong(em("Tsplit_max."))),
                   h3("For each simulation in the ", strong("SC"), " and ", strong("AM"), " models, the time of secondary contact between lines", strong(em("(Tsc),")), " or old migration stop times", strong(em("(Tam)")), "are drawn uniformly between", strong(em("Tsplit_min")), "and", strong(em("Tsplit_sampled.")))
               ),
               
              # box(title = h2("Migration rates"), icon = NULL, width = NULL, solidHeader = TRUE, status = "success", height = 250,
              #boxPlus(title = h2("Migration rates") %>% helper(type = "inline", title = "Inline Help", content = c("This helpfile is defined entirely in the UI!", "This is on a new line.", "This is some <b>HTML</b>."), size = "s"), height = 250, width = NULL, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
              boxPlus(title = h2("Migration rates"), height = 250, width = NULL, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
                      
                                      fluidRow(
                       column(width=5, numericInput("M_min", label = h5('min'), value = 0.4)),
                       column(width=5, numericInput("M_max", label = h5('max'), value = 20))
                       )
              ),
              
              boxPlus(
                  title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = NULL,
                  boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                  enable_label = TRUE, label_text = "Informations about migration rates", label_status = "success",
                  
                  h3("Migration rates are expressed in", strong(em("4.Ne.m,")), "where", strong(em("m")), "is the fraction of each subpopulation made up of new migrants each generation.")
              )
        )
    ),
    
    fluidRow(
        box("", width = 12, solidHeader = TRUE, status = "info",
            prettyCheckbox(inputId = "check_prior", shape = "round", value = FALSE,
                           label = strong("Please check/valid your choices"), icon = icon("check"),
                           animation = "tada", status = "success", bigger = TRUE)
        )
    )
)



run_ABC <- dashboardBody(
    # PRINT INFOX BOXES OF CHECKING
    fluidRow(
            boxPlus(
#                shiny::tags$h3("Checked options"),
                width = 12,
                uiOutput('check_upload_info'),
                uiOutput('check_filtering_info'),
                uiOutput('check_populations_info'),
                uiOutput('check_prior_info')
            )
        ),
    
    # PRINT INPUT
    fluidRow(
        column(width = 12,
            boxPlus(
                title = h2("Information summary"), width = NULL, icon = "fa fa-heart", solidHeader = TRUE, gradientColor = "teal",
                boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                enable_label = TRUE, label_text = "Please check the following information", label_status = "success",

                tableOutput("parameters")
            )
        ),
    
    # RUN ABC
    uiOutput("run_ABC")
    )
)

upload_results <- dashboardBody(
  # upload results
  fluidRow(NULL, soldHeader = TRUE, status ="danger",
      boxPlus(title = h2("Results to upload (i.e, fastABC's archived output)"), height = 200,  width = 12, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
      fileInput("results", label = NULL)
    )
  )
)


distribution_statistics <- dashboardBody(
#  boxPlus(title = h2("Distribution of observed statistics"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
          
#    fluidRow(
#      # Sidebar panel for inputs
#      column(width=3,
 #     # Input: Slider for the number of bins
  #    #sliderInput(inputId = "nObsStats", label = h3("Number of observed statistics to plot:"), min = 1, max = 13, value = 1, step=1)
   #   selectInput("obs_stat_to_plot", label = h3("Statistics to plot"), choices = list("Summary of the SFS" = 1, "Diversity" = 2, "Tajima's D" = 3, "Differentiation and divergence" = 4), selected = 1)
    #  ),
            
#    # Main panel for displaying outputs
 #   column(width=9,
  #  # Output: density

   # plotlyOutput("plot_obs_stats")
    #)
 # )
#))
  tags$head(
    tags$style(type='text/css', 
               ".nav-tabs {font-size: 18px} ")),    
  
  tabsetPanel(type = "tabs",
              tabPanel("Summarized jSFS", plotlyOutput("plot_obs_stats_sites")),
              tabPanel("Polymorphism", plotlyOutput("plot_obs_stats_diversity")),
              tabPanel("Tajima's D", plotlyOutput("plot_obs_stats_tajima")),
              tabPanel("Differentiation and divergence", plotlyOutput("plot_obs_stats_divergence"))
  )
)


continuum_divergence <- dashboardBody(
  tags$head(tags$script('
                        var dimension = [0, 0];
                        $(document).on("shiny:connected", function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                        $(window).resize(function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                        ')),
  
  #oxPlus(title = h2("Speciation along a continuum of divergence"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
  #plotlyOutput("plot_greyzone")
  # align="middle" height="auto" width="100%" margin="0 auto
  htmltools::div(style = "display:inline-block", plotlyOutput("plot_greyzone", width = "auto"))
  
  
  
)


informations <- dashboardBody(
    fluidRow(
      column(width = 12,
#        box(title = h2("Citations"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
        boxPlus(title = h2("Citations"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
            h3("Please, in case of publication of a study using fastABC, do not forget to quote the following references:"),
            h3(code('Csilléry, Katalin, Olivier François, and Michael GB Blum. "abc: an R package for approximate Bayesian computation (ABC)." Methods in ecology and evolution 3.3 (2012): 475-479.')),
            h3(code('Pudlo, Pierre, Jean-Michel Marin, Arnaud Estoup, Jean-Marie Cornuet, Mathieu Gautier, and Christian P. Robert. "Reliable ABC model choice via random forests." Bioinformatics 32, no. 6 (2015): 859-866.')),
            h3(code('Roux, Camille, Christelle Fraisse, Jonathan Romiguier, Yoann Anciaux, Nicolas Galtier, and Nicolas Bierne. "Shedding light on the grey zone of speciation along a continuum of genomic divergence." PLoS biology 14, no. 12 (2016): e2000234.'))
        ),
        
#        box(title = h2("Acknowledgment"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
        boxPlus(title = h2("Acknowledgment"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
            h3("Please, if you use this online version of fastABC, do not forget to recognize and acknowledge the free provision of calculation cores by France Bioinformatique"),
            h3(code("The demographic inferences were conducted on the IFB Core Cluster which is part of the National Network of Compute Resources (NNCR) of the", a(span(strong("Institut Français de Bioinformatique (IFB)."), style = "color:teal"), href="https://www.france-bioinformatique.fr/fr", target="_blank")))
        ),
        
#        box(title = h2("Partners"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
        boxPlus(title = h2("Partners"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
            mainPanel(htmlOutput("logos"))
        )
      )
    )
)

ui <- dashboardPage(
    dashboardHeader(title = "menu fastABC"
      ),
    dashboardSidebar(
      tags$head(
        tags$style(HTML(".main-sidebar { font-size: 40px; }")) #change the font size to 20
      ),
      
      sidebarMenu(
#            style = "position: fixed; overflow: visible;",
            menuItem(("Welcome"), tabName = "welcome", icon = icon("door-open")),
            menuItem(('ABC'), tabName = "ABC", icon = icon("industry"),
              menuSubItem(("Upload data"), tabName = "upload", icon = icon("cloud-upload")),
              menuSubItem(("Data filtering"), tabName = "filtering", icon = icon("bath")),
              menuSubItem(("Populations/species"), tabName = "populations", icon = icon("users-cog")),
              menuSubItem(("Prior distributions"), tabName = "simulations", icon = icon("dice")),
              menuSubItem(("Run ABC"), tabName = "run_abc", icon = icon("microchip"))
            ),
            menuItem(("Results visualization"), tabName = "Results_visualization", icon = icon("chart-pie"),
                     menuSubItem(("Upload results"), tabName = "upload_results", icon = icon("cloud-upload")),
                     menuSubItem(("Distribution of statistics"), tabName = "distribution_statistics", icon = icon("chart-bar")),
                     menuSubItem(("Continuum of divergence"), tabName = "continuum_divergence", icon = icon("arrow-alt-circle-right"))
            ),
            menuSubItem(("Informations"), tabName = "information", icon = icon("info-circle"))
        )
    ),
    
    dashboardBody(
      tags$head(tags$style(
        HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden;}')
      )),
      
      setShadow(class = "box"),
      shinyDashboardThemes(
            theme = "poor_mans_flatly"
        ),
        
        tags$head( 
            tags$style(HTML(".main-sidebar { font-size: 40px; }")) #change the font size to 20
        ),
        
        tabItems(
            # Welcome
            tabItem(tabName = "welcome",
                welcome_page
            ),
            
            # Upload data
            tabItem(tabName = "upload",
                upload_data
            ),
        
            #  Filtering
            tabItem(tabName = "filtering",
                filtering
            ),
            
            # Populations
            tabItem(tabName = "populations",
                populations
            ),
 
            # Simulations
            tabItem(tabName = "simulations",
                prior
            ),

            # Run the BAC inferences
            tabItem(tabName = "run_abc",
                run_ABC
            ),
            
            # Upload the fastABC's results
            tabItem(tabName = "upload_results",
                upload_results
            ),
            
            # Plot the distributions of statistics
            tabItem(tabName = "distribution_statistics",
                distribution_statistics
            ),
            
            # Plot the distributions of statistics
            tabItem(tabName = "continuum_divergence",
                    continuum_divergence
            ),
            
            # Informations
            tabItem(tabName = "information",
                informations
            )
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output, session = session) {
    options(shiny.maxRequestSize=4000*1024^2)
    #  WELCOME
    output$welcome_picture <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/dag_2pops.pdf.png align="middle" height="auto" width="100%" margin="0 auto">'
        )
      }
    )
    
    output$models_picture <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/models.png align="middle" height="auto" width="100%" margin="0 auto">'
          )
        }
      )
    
    output$model_comparisons <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/model_comparisons.png align="middle" height="auto" width="100%" margin="0 auto">'
        )
      }
    )
    
    output$homo_hetero <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/homo_hetero.png align="middle" height="auto" width="100%" margin="0 auto">'
        )
      }
    )
    
    
    ## example of genomic heterogeneity
    output$genomic_hetero <- renderPlot({
      par(las=1)
      y_points = dbeta(0:100/100, input$alpha, input$beta)
      x_points = 0:100/100 * input$Ne
      plot(x_points, y_points, type='l', xlab = expression(paste("Genomic distribution of ", italic('Ne'), sep=" ")), ylab='density', main=expression(italic("Example of genomic distributions that fastABC will try to infer")), col="white", cex.main = 1.5, cex.axis = 1.5, cex.lab=1.5, xlim=c(min(c(x_points, input$Ne*1.2)), max(c(x_points, input$Ne*1.2))))
      
      x_points = c(0, x_points, max(x_points), 0)
      y_points = c(0, y_points, 0, 0)
      polygon(x_points, y_points, border = 'NA', col="#b2df8a")
      abline(v=input$Ne, lwd=4, col="#fee08b")
    })
    
    
    #  UPLOAD DATA
    ## GET THE SUMMARY STATS ABOUT THE UPLOADED FILE
        ## list of species
        list_species = reactive({
            if(is.null(input$infile)){return ()}
          withProgress(message = 'Getting the species', detail = NULL, value = 0, {
            incProgress(1/2)
            return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f2 | sort -u", sep=" "), intern = T))
            incProgress(1/2)
          })
        })
        output$list_species <- renderDataTable(data.frame("species" = list_species()))

        ## list of individuals
        list_individuals = reactive({
            if(is.null(input$infile)){return ()}
          withProgress(message = 'Getting the individuals', detail = NULL, value = 0, {
            incProgress(1/2)
            return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f3 | sort -u", sep=" "), intern = T))
            incProgress(1/2)
          })
        })
        output$list_individuals <- renderDataTable(data.frame("individuals" = list_individuals()))

        ## list of loci
        list_loci = reactive({
            if(is.null(input$infile)){return ()}
          withProgress(message = 'Getting the loci', detail = NULL, value = 0, {
            incProgress(1/2)
            return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f1 | sort -u | cut -d'>' -f2", sep=" "), intern = T))
            incProgress(1/2)
          })
        })
        output$list_loci <- renderDataTable(data.frame("loci" = list_loci()))

            ## Summary Stats code ##
        # this reactive output contains the summary of the dataset and display the summary in table format
        nSpecies = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f2 | sort -u", sep=" "), intern = T))})
        nIndividuals = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f3 | sort -u", sep=" "), intern = T))})
        nLoci = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f1 | sort -u", sep=" "), intern = T))})
        output$general_informations <- renderTable({
            if(is.null(input$infile)){return ()}
          withProgress(message = 'Producing the table of informations', detail = NULL, value = 0, {
            incProgress(1/2)
            return(data.frame("nSpecies" = length(list_species()), "nIndividuals" = length(list_individuals()), "nLoci" = length(list_loci())))
            incProgress(1/2)
          })
        })

        ## MainPanel tabset renderUI code ##
        # the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
        # Until the file is loaded, app will not show the tabset.
        output$upload <- renderUI({
            if(is.null(input$infile)) {return(loadingState())}
            else
                tabsetPanel(
                    tabPanel("General informations", tableOutput("general_informations")),
                    tabPanel("List of individuals", dataTableOutput("list_individuals")),
                    tabPanel("List of populations or species", dataTableOutput("list_species")),
                    tabPanel("List of loci", dataTableOutput("list_loci"))
                )
        })

        ## Number of ABC analysis to perform
      #  observeEvent(input$number_of_ABC_validation, {
        observe(if(input$check_upload){
            shinyjs::disable("number_of_ABC")
        })
        
        
    # POPULATIONS/SPECIES
        output$input_names_ui <- renderUI({
            if(is.null(input$infile)) {return()}
            nspecies = as.integer(input$nspecies)
            lapply(1:nspecies, function(i){
                selectInput(paste0("name", LETTERS[i]), label = paste0("name of species/population ", LETTERS[i]), choices = list_species(), selected = paste0("name", LETTERS[i]))
            })
        })
        
        output$input_names_outgroup_ui <- renderUI({
            if(is.null(input$infile)) {return()}
            presence_outgroup = input$presence_outgroup
            if(presence_outgroup == 'no'){
                # nameOutgroup = 'NA' # CONFIG_YAML
            }else{
                selectInput("nameOutgroup", label = "name of the outgroup species", choices = list_species())
                # nameOutgroup = 'name specified by the selectInput # CONFIG_YAML
            }
        })
        
    # PRINT INPUT
    output$parameters <- renderTable({
        if(is.null(input$infile)) {return()}
        mail_address = input$mail_address # email address
        config_yaml = "not to be specified" # name of the config_yaml used by snakemake
        infile = input$infile$name # name of the input fasta file
        region = input$region # coding or noncoding
        nspecies = input$nspecies # number of species to simulate
        if(nspecies == 1){
            species_names = c(input$nameA)
            species_names_row = c('nameA')
        }else{
            if(nspecies == 2){
                species_names = c(input$nameA, input$nameB)
                species_names_row = c('nameA', "nameB")
            } else{
                if(nspecies == 4){
                    species_names = c(input$nameA, input$nameB, input$nameC, input$nameD) # name of the simulated species
                    species_names_row = c('nameA', "nameB", "nameC", "nameD")
                }else{
                    species_names = "NA"
                    species_names_row = "NA"
                }
            }
        }
        
        if(input$presence_outgroup == 'yes'){
            nameOutgroup = input$nameOutgroup
        }else{
            nameOutgroup = "NA"
        }
        
        Lmin = input$Lmin
        nMin = input$nMin
        mu = input$mu
        rho_over_theta = input$rho_over_theta
        N_min = input$N_min
        N_max = input$N_max
        Tsplit_min = input$Tsplit_min
        Tsplit_max = input$Tsplit_max
        M_min = input$M_min
        M_max = input$M_max
        
        res = matrix(c(mail_address, config_yaml, infile, region, nspecies, species_names, nameOutgroup, Lmin, nMin, mu, rho_over_theta, N_min, N_max, Tsplit_min, Tsplit_max, M_min, M_max), ncol = 1)
        
        row.names(res) = c("user's email address", "config_yaml", "infile", "region", "nspecies", species_names_row, "nameOutgroup", "Lmin", "nMin", "mu", "rho_over_theta", "N_min", "N_max", "Tsplit_min", "Tsplit_max", "M_min", "M_max")
        colnames(res) = c("entries")
        
        res        
    }, rownames = TRUE, colnames = TRUE)
    
    # RUN ABC
    ## only show the action button RUN ABC if a file is uploaded and 4 checkings were made
    output$run_ABC <- renderUI({
    if(is.null(input$infile)==FALSE && input$check_upload == TRUE && input$check_filtering == TRUE && input$check_populations == TRUE && input$check_prior == TRUE){
        a <-column(width = 12,
               boxPlus(
                   title = h2("Run ABC"), width = NULL, icon = "fa fa-heart", solidHeader = TRUE, gradientColor = "teal",
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "Are you ready?", label_status = "danger",
                   h3("Submission of the ABC workflow"),
                   actionButton("runABC", label = "Run the ABC", size = 'md', width = '100%', fullwidth = TRUE),
                   h3("Number of submitted analysis (you can change the studied populations/species and the priors between 2 analysis):"),
                   verbatimTextOutput('nClicks'),
                   hr(),
                   h3("Timestamp of the last submitted analysis:"),
                   verbatimTextOutput('time_stamp'),
                   hr(),
                   h3("Snakemake command line:"),
                   verbatimTextOutput("snakemake_command")
               )
            )
        }else{return()}
    })
    
    ## print the number of clicks on the "run the ABC" button
    output$nClicks <- renderText({ input$runABC })
    
    ## removing the "Run the ABC" button after clicking on it
    observeEvent(input$runABC, {if (input$runABC == input$number_of_ABC)  removeUI(selector='#run_ABC', immediate=TRUE)}, autoDestroy=TRUE)
    
    ## get a time stamp when clicking on the "Run the ABC" button
    time_stamp <- reactiveVal(0)
    observeEvent( input$runABC, {time_stamp(system('echo $(mktemp -d -t XXXXXXXXXX | cut -d"/" -f3)', intern=T))})
    output$time_stamp <- renderText({time_stamp()})

    ## snakemake command
    snakemake_command <- reactiveVal(0)
    observeEvent( input$runABC, {snakemake_command(paste('snakemake -p -j 999 --snakefile ../2pops/Snakefile --configfile config_', time_stamp(), '.yaml --cluster-config ../cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}"', sep=''))})
    output$snakemake_command <- renderText({snakemake_command()})
        
    ## Check upload
    output$check_upload_info <- renderUI({
        if(input$check_upload == FALSE) {
            a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_upload == TRUE){
            a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check filtering
    output$check_filtering_info <- renderUI({
        if(input$check_filtering == FALSE) {
            a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("bath"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_filtering == TRUE){
            a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("bath"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check populations
    output$check_populations_info <- renderUI({
        if(input$check_populations == FALSE) {
            a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_populations == TRUE){
            a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check prior
    output$check_prior_info <- renderUI({
        if(input$check_prior == FALSE) {
            a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("dice"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_prior == TRUE){
            a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("dice"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## RESULT VISUALIZATION
    locus_spe <- reactive({
      fileName = input$results
      
      if (is.null(fileName))
        return(NULL)
      
      untar(fileName$datapath, exdir = getwd())

      rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
      
      locus_spe_name = paste(rootName, "/modelComp/locus_specific_modelComp.txt", sep='')
      print(locus_spe_name)
      #read the table
      locus_spe = read.table(locus_spe_name, h=T)
      
      # delete the untar results
      system(paste('rm -rf ', rootName, sep=''))
      
      # return the read object
      return(locus_spe)
    })
    
    # GRAPH SITES
    output$plot_obs_stats_sites <- renderPlotly({
      # number of loci
      nLoci = nrow(locus_spe())
      
      f <- list(
        family = "Arial",
        size = 20
      )
      axis_x <- list(
        title = "",
        titlefont = f,
        tickfont = list(size = 20)
      )     
      axis_y <- list(
        title = "Proportion of sites",
        titlefont = f,
        tickfont = list(size = 20)
      )
      
      statistics_obs_sites = c(locus_spe()$sf_avg, locus_spe()$sxA_avg, locus_spe()$sxB_avg, locus_spe()$ss_avg)
      statistics_names_sites = rep(c("proportion of\nfixed differences\nbetween A and B", "proportion of\npolymorphic sites\nexclusive to A", "proportion of\npolymorphic sites\nexclusive to B", "proportion of\nshared polymorphic sites\nbetween A and B"), each = nLoci)
      data_obs_sites = data.frame(statistics_obs_sites, statistics_names_sites)
      head(data_obs_sites)
      graph_sites = plot_ly(data_obs_sites, y=~statistics_obs_sites, x=~statistics_names_sites, color=~statistics_names_sites, type="violin", box = list( visible = T ), width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)))     
      return(graph_sites)
    })
    
    # GRAPH DIVERSITY
    output$plot_obs_stats_diversity <- renderPlotly({
      # number of loci
      nLoci = nrow(locus_spe())
      
      f <- list(
        family = "Arial",
        size = 20
      )
      axis_x <- list(
        title = "",
        titlefont = f,
        tickfont = list(size = 20)
      )
      axis_y <- list(
        title = "Index of diversity per site",
        titlefont = f,
        tickfont = list(size = 20)
      )
      statistics_obs_diversity = c(locus_spe()$piA_avg, locus_spe()$piB_avg, locus_spe()$thetaA_avg, locus_spe()$thetaB_avg)
      statistics_names_diversity = rep(c("diversity in sp. A\n(measured by pi)", "diversity in sp. B\n(measured by pi)", "diversity in sp. A\n(measured by Watterson's theta)", "diversity in sp. B\n(measured by Watterson's theta)"), each = nLoci)
      data_obs_diversity = data.frame(statistics_obs_diversity, statistics_names_diversity)
      
      graph_diversity = plot_ly(data_obs_diversity, y=~statistics_obs_diversity, x=~statistics_names_diversity, color=~statistics_names_diversity, type="violin", box = list( visible = T ), width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)))
      return(graph_diversity)
    })
    
    
    # GRAPH TAJIMA
    output$plot_obs_stats_tajima <- renderPlotly({
      # number of loci
      nLoci = nrow(locus_spe())
      f <- list(
        family = "Arial",
        size = 20
      )
      axis_x <- list(
        title = "",
        titlefont = f,
        tickfont = list(size = 20)
      )
      axis_y <- list(
        title = "Tajima's D",
        titlefont = f,
        tickfont = list(size = 20)
      )
      statistics_obs_tajima = c(locus_spe()$DtajA_avg, locus_spe()$DtajB_avg)
      statistics_names_tajima = rep(c("Tajima's D in sp. A", "Tajima's D in sp. B"), each = nLoci)
      data_obs_tajima = data.frame(statistics_obs_tajima, statistics_names_tajima)
      
      graph_tajima = plot_ly(data_obs_tajima, y=~statistics_obs_tajima, x=~statistics_names_tajima, color=~statistics_names_tajima, type="violin", box = list( visible = T ), width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(2)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)))
      return(graph_tajima)
    })
    
    
    # GRAPH DIVERGENCE
    output$plot_obs_stats_divergence <- renderPlotly({
      # number of loci
      nLoci = nrow(locus_spe())
      f <- list(
        family = "Arial",
        size = 30
      )
      axis_x <- list(
        title = "",
        tickfont = list(size = 20)
      )
      
      axis_y <- list(
        title = "Measure of divergence/differentiation",
        titlefont = f,
        tickfont = list(size = 20)
      )
      
      statistics_obs_divergence = c(locus_spe()$divAB_avg, locus_spe()$netdivAB_avg, locus_spe()$FST_avg)
      statistics_names_divergence = rep(c("raw divergence between\nsp. A and sp. B", "net divergence between\nsp. A and sp. B", "Fst between\nsp. A and sp. B"), each = nLoci)
      data_obs_divergence = data.frame(statistics_obs_divergence, statistics_names_divergence)
      
      graph_divergence = plot_ly(data_obs_divergence, y=~statistics_obs_divergence, x=~statistics_names_divergence, color=~statistics_names_divergence, type="violin", box = list( visible = T ), width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(3)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)) )
      return(graph_divergence)
    })
    
    
    output$plot_greyzone <- renderPlotly({
      x = read.table("popPhyl.txt", h=T)
      col = c(grey(0.25), 'turquoise', 'purple', 'red')
      #col = c(grey(0.75), rev(viridis(5))[3:5])
      pmig_HH = x$Pongoing_migration_Mhetero_Nhetero 
      proba_migration = pmig_HH
      seuil1 = 0.6419199
      seuil2 = 0.1304469
      
      model = rep('ambiguous', nrow(x))
      model[which(x$Pongoing_migration_Mhetero_Nhetero>=seuil1)] = "migration"
      model[which(x$Pongoing_migration_Mhetero_Nhetero<seuil2)] = "isolation"
      
      divergence = log10(x$netdivAB_avg)
      
      piA = round(x$piA_avg, 5)
      piB = round(x$piB_avg, 5)
      
      pattern=c("Mhetero_Nhetero", "Hetero")
      selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(x))))
      status = rep('ambiguous', nrow(x))
      heteroM = apply(x[, selectedCol], FUN="sum", MARGIN=1)
      status[which(pmig_HH>= seuil1 & heteroM >= seuil1)] = "semi-isolated species"
      status[which(pmig_HH>= seuil1 & heteroM < seuil1)] = "populations"
      status[which(pmig_HH<= seuil2)] = "species"
      
      species_A = x$spA
      species_B = x$spB
      
      author = rep('camille.roux@univ-lille.fr', length(species_A))
      res = data.frame(divergence, model, status, proba_migration, species_A, species_B, piA, piB, author)
      
      f=list(
        family = "Arial",
        size = 26,
        color = "black"
      )
      
      f2=list(
        family = "Arial",
        size = 20,
        color = "black"
      )
      
      f_legend=list(
        family = "Arial",
        size = 20,
        color = "black",
        color = "#000"
      )
      
      xlab = list(title='divergence (log10)',
                  titlefont=f,
                  tickfont=f2,
                  tickvals=c(0, -1, -2, -3, -4, -5),
                  #ticktext=c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5")))
                  ticktext=c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
      )
      
      ylab = list(title='probability of ongoing migration',
                  titlefont=f,
                  tickfont=f2
      )
      
      #p=plot_ly(data=res, x=~divergence, y=~proba_migration, color=~status, colors=col, marker=list(size=20), text = ~paste("species A: ", species_A, '<br>species B: ', species_B, "<br>net neutral divergence: ", round(10**divergence, 5), "<br>Probability of migration: ", round(proba_migration, 3), '<br>piA: ', piA, '<br>piB: ', piB, '<br><br>author: ', author), hoverinfo='text') %>% layout(xaxis=xlab, yaxis=ylab, legend=list(orientation = 'h', y=1.05, font=f_legend))
      p=plot_ly(data=res, x=~divergence, y=~proba_migration, color=~status, colors=col, marker=list(size=20), text = ~paste("species A: ", species_A, '<br>species B: ', species_B, "<br>net neutral divergence: ", round(10**divergence, 5), "<br>Probability of migration: ", round(proba_migration, 3), '<br>piA: ', piA, '<br>piB: ', piB, '<br><br>author: ', author), hoverinfo='text', width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2])) %>% layout(xaxis=xlab, yaxis=ylab, legend=list(orientation = 'h', y=1.05, font=f_legend))
      #htmlwidgets::saveWidget(p, "figure_greyzone.html") # HTML
      #webshot::webshot("figure_greyzone.html", "figure_greyzone.pdf") # PDF -> Margin problem (cut!)
      return(p)
      
    })
      
    ## INFORMATIONS
    ### LOGOS
    output$logos <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/logos.png align="middle" height="auto" width="100%" margin="0 auto">'
        )
      }
      )
}

# Run the application 
shinyApp(ui = ui, server = server)

