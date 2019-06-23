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
library(shinymaterial)


# welcome
welcome_page <- dashboardBody(
  fluidRow(
    box(
      title = h2("Overview"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
        h3(strong("fastABC"), "is a DNA sequence analysis workflow to study the demographic history of sampled populations or species by using Approximate Bayesian Computations."),
        h3("From a single uploaded input file containing sequenced genes or DNA fragments,", strong("fastABC"), "will:"),
        h3(strong("1."), "simulate different models/scenarios."),
        h3(strong("2."), "select the best model using an ABC approach based on", a(span(strong("random forests."), style = "color:teal"), href="https://cran.r-project.org/web/packages/abcrf/index.html", target="_blank")),
        h3(strong("3."), "estimate the parameters of the best model using a", a(span(strong("neural network"), style = "color:teal"), href="https://cran.r-project.org/web/packages/abc/index.html", target="_blank"), "approach."),
        h3(strong("4."), "measure the robustness of the analyses.", strong("fastABC"), "is transparent on the ability of its inferences to reproduce the observed data."),
        hr(),
      
        h3("FastABC is composed of two elements:"),
        h3(strong("1."), "A web interface developed in", a(span(strong("Shiny,"), style = "color:teal"), href="https://shiny.rstudio.com/", target="_blank"), "which will execute"),
        h3(strong("2."), "a workflow managed by", a(span(strong("Snakemake."), style = "color:teal"), href="https://snakemake.readthedocs.io/en/stable/", target="_blank")),
        h3("The code is fully open-source, freely distributed on", a(span(strong("GitHub"), style = "color:teal"), href="https://github.com/popgenomics/ABConline", target="_blank"), "and can be immediately redeployed on any cluster using", a(span(strong("SLURM"), style = "color:teal"), href="https://slurm.schedmd.com/documentation.html", target="_blank"), "thanks to a", a(span(strong("Singularity"), style = "color:teal"), href="https://sylabs.io/docs/#doc-3.2", target="_blank"), "image."),
        
        h3("This redeployment allows the user to modify the models to be compared, to add summary statistics, etc..."),
        h3("The workflow can be simply executed from the command line without going through the web interface."),
        
        h3("However, if desired, the web interface can also be freely hosted and linked to any cluster.")
      )
    ),
  
  fluidRow(
    box(
      title = h2("Compared models"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
                mainPanel(htmlOutput("models_picture"))
      
    ),
    
    box(
      title = h2("Model comparisons for 2 populations/species"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
      mainPanel(htmlOutput("model_comparisons"))
    ),
    
    box(
      title = h2("Snakemake pipeline"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
                mainPanel(htmlOutput("welcome_picture"))
    )
  )
)

# upload
upload_data <- dashboardBody(
        tags$head(tags$style(HTML("a {color: black}"))),
        fluidRow(
            box(
                title = h2("Number of ABC analysis to run"), height = 250, width = 6, solidHeader = TRUE, background = NULL, status = "primary",
                shinyjs::useShinyjs(),
                #numericInput(inputId = "number_of_ABC", label = NULL,  value = 1, min = 1, max = 5, step = 1),
                selectInput("number_of_ABC", label = h4("1 to 5 ABC analyses can be performed from the same input file"), choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5), selected = 1)
                
               # actionButton("number_of_ABC_validation", "Validate")
            ),
            
            box(
                title = h2("Email address"), height = 250,  width = 6, solidHeader = TRUE, status = "primary",
                textInput("mail_address", label = h4("address to receive the download link of the results"), value = "user@gmail.com")
            )
        ),

     fluidRow(NULL, soldHeader = TRUE, status ="danger",
              box(
                  title = h2("Sequence Alignment Upload"), height = 200,  width = 6, solidHeader = TRUE, status = "success",
                  fileInput("infile", label = NULL)
              ),
              
             
             box(
                 title = h2("Genomic regions"), height = 200,  width = 6, solidHeader = TRUE, status = "warning",
                 #                radioButtons("region", label = NULL, selected = "noncoding", choices = list("coding" = "coding", "non coding" = "noncoding"))
                 prettyRadioButtons("region", label = NULL, shape = "round", status = "warning", fill = TRUE, inline = TRUE, animation = "pulse", bigger = TRUE,
                                    selected = "noncoding", choices = list("coding" = "coding", "non coding" = "noncoding"))
             )
        ),     
        
        fluidRow(align="left",
            box(
                title = h2("Input file"), width = 6, solidHeader = TRUE, background = "green", status = "success",
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
            
            box(
                title = h2("Informations about the uploaded file"),  width = 6, solidHeader = TRUE, status = "success",
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
               box(
                   title = h2("Maximum proportion of N"), width = NULL, solidHeader = TRUE, status = "primary", height = 225,
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
               box(
                   title = h2("Minimum gene length"), width = NULL, solidHeader = TRUE, status = "success", height = 225,
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
                    h3("Lmin = (number of ", strong("monomorphic positions"), "that can be considered) + (number of ", strong("biallelic positions"), "that can be considered)."),
                    hr(),
                    h3(a(span(strong("Example", style = "color:green")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/Lmin.png", target="_blank"))
               )
        ),
        
        column(width = 4,
               box(
                   title = h2("Minimum number of sequences"), width = NULL, solidHeader = TRUE, status = "warning", height = 225,
                   numericInput("nMin", label = NULL, value = 12)
               ),
               boxPlus(
                   title = h3("nMin"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "EXPLANATIONS", label_status = "warning",
                   h3(strong("fastABC"), " starts for each gene by eliminating individual sequences containing too many N ", span(strong("(max_N_tolerated; blue box)", style = "color:blue")), "."),
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
        box(
            title = h2("Number of populations/species"), width = 4, solidHeader = TRUE, status = "primary", height = 600,
#            radioButtons("nspecies", label = h3("Number of gene pools"), choices = list("One gene pool" = 1, "Two gene pools" = 2, "Four gene pools" = 4), selected = 2),
            prettyRadioButtons("nspecies", label = h3("Number of gene pools"), shape = "round", status = "primary", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
                               choices = list("One gene pool" = 1, "Two gene pools" = 2, "Four gene pools" = 4), selected = 2),
            uiOutput("input_names_ui")
        ),
        
        box(
            title = h2("Outgroup species"), width = 4, solidHeader = TRUE, status = "warning", height = 600,
#            radioButtons("presence_outgroup", label = h3("Presence of an outgroup"), choices = list("no" = "no", "yes" = "yes"), selected = "no"),
            prettyRadioButtons("presence_outgroup", label = h3("Presence of an outgroup"), shape = "round", status = "warning", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
                               choices = list("no" = "no", "yes" = "yes"), selected = "no"),
            uiOutput("input_names_outgroup_ui")
        )
    ),

    fluidRow(
        box("", width = 8, solidHeader = TRUE, status = "info",
            prettyCheckbox(inputId = "check_populations", shape = "round", value = FALSE,
                           label = strong("Please check/valid your choices"), icon = icon("check"),
                           animation = "tada", status = "success", bigger = TRUE)
        )
    )
)


prior <- dashboardBody(
    fluidRow(
        column(width = 4,
               box(
               #    title = strong("Mutation and recombination"), icon = "fas fa-bolt", width = 12, solidHeader = TRUE, gradientColor = "olive", boxToolSize = "lg", footer_padding = TRUE, collapsible = FALSE, collapsed = FALSE, closable = FALSE,
                   title = h2("Mutation and recombination"), width = NULL, solidHeader = TRUE, status = "primary", height = 250,
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
               
               box(
                    #title = strong("Ne (# of diploid individuals)"), icon = "fas fa-users", width = 12, solidHeader = TRUE, gradientColor = "purple", boxToolSize = "lg", footer_padding = TRUE, collapsible = FALSE, collapsed = FALSE, closable = FALSE,
                   
                    title = h2("Population sizes"), width = NULL, solidHeader = TRUE, status = "danger", height = 250,
                    fluidRow(
                        column(width=5, numericInput("N_min", label = h5('min'), value = 100)),
                        column(width=5, numericInput("N_max", label = h5('max'), value = 1000000))
                   )
               ),
               
               boxPlus(
                   title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "danger",
                   boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                   enable_label = TRUE, label_text = "Informations about population sizes", label_status = "danger",
                   
                   h3("The effective population sizes ", em(strong("Ne")), "is the number of diploid individuals within current and ancestral species/populations."),
                   h3("In the", strong("ABC"), "simulations,", em(strong("Ne")), "is independent between all current and ancestral species/populations.")
               )
        ),
        
        column(width = 4,
               box(
                   #title = strong("Time of split (# of generations)"), icon = "fas fa-history", width = 12, solidHeader = TRUE, gradientColor = "teal", boxToolSize = "lg", footer_padding = TRUE, collapsible = FALSE, collapsed = FALSE, closable = FALSE,
                   title = h2("Time of split"), width = NULL, solidHeader = TRUE, status = "warning", height = 250,
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
               
               box(
                   #title = strong("Migration rate (4.N.m per generation)"), icon = "fas fa-helicopter", width = 12, solidHeader = TRUE, gradientColor = "teal", boxToolSize = "lg", footer_padding = TRUE, collapsible = FALSE, collapsed = FALSE, closable = FALSE,
                   
                   title = h2("Migration rates"), icon = NULL, width = NULL, solidHeader = TRUE, status = "success", height = 250,
                   fluidRow(
                       column(width=5, numericInput("M_min", label = h5('min'), value = 0.4)),
                       column(width=5, numericInput("M_max", label = h5('max'), value = 20))
                       )
              ),
              
              boxPlus(
                  title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = NULL,
                  boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
                  enable_label = TRUE, label_text = "Informations about migration rates", label_status = "success",
                  
                  h3("Migration rates are expressed in", strong(em("4.N.m,")), "where", strong(em("m")), "is the fraction of each subpopulation made up of new migrants each generation.")
              )
        )
    ),
    
    fluidRow(
        box("", width = 8, solidHeader = TRUE, status = "info",
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
                width = 10,
                uiOutput('check_upload_info'),
                uiOutput('check_filtering_info'),
                uiOutput('check_populations_info'),
                uiOutput('check_prior_info')
            )
        ),
    
    # PRINT INPUT
    fluidRow(
        column(width = 10,
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


informations <- dashboardBody(
    fluidRow(
      column(width = 10,
        box(title = h2("Citations"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
            h3("Please, in case of publication of a study based on fastABC, do not forget to quote the following references:"),
            h3(code('Csilléry, Katalin, Olivier François, and Michael GB Blum. "abc: an R package for approximate Bayesian computation (ABC)." Methods in ecology and evolution 3.3 (2012): 475-479.')),
            h3(code('Pudlo, Pierre, Jean-Michel Marin, Arnaud Estoup, Jean-Marie Cornuet, Mathieu Gautier, and Christian P. Robert. "Reliable ABC model choice via random forests." Bioinformatics 32, no. 6 (2015): 859-866.')),
            h3(code('Roux, Camille, Christelle Fraisse, Jonathan Romiguier, Yoann Anciaux, Nicolas Galtier, and Nicolas Bierne. "Shedding light on the grey zone of speciation along a continuum of genomic divergence." PLoS biology 14, no. 12 (2016): e2000234.'))
        ),
        
        box(title = h2("Acknowledgment"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
            h3("Please, if you use this online version of fastABC, do not forget to recognize the free provision of calculation cores by France Bioinformatique"),
            h3(code("The demographic inferences were conducted on the IFB Core Cluster which is part of the National Network of Compute Resources (NNCR) of the", a(span(strong("Institut Français de Bioinformatique (IFB)."), style = "color:teal"), href="https://www.france-bioinformatique.fr/fr", target="_blank")))
        )
      )
    )
)

ui <- dashboardPage(
    dashboardHeader(title = "fastABC"),
    dashboardSidebar(
        sidebarMenu(
            style = "position: fixed; overflow: visible;",
            menuItem(h3("Welcome"), tabName = "welcome", icon = icon("door-open", "fa-2x")),
            menuItem(h3("Upload data"), tabName = "upload", icon = icon("cloud-upload", "fa-2x")),
            menuItem(h3("Data filtering"), tabName = "filtering", icon = icon("bath", "fa-2x")),
            menuItem(h3("Populations/species"), tabName = "populations", icon = icon("users-cog", "fa-2x")),
            menuItem(h3("Prior distributions"), tabName = "simulations", icon = icon("dice", "fa-2x")),
            menuItem(h3("Run ABC"), tabName = "run_abc", icon = icon("microchip", "fa-2x")),
            menuItem(h3("Informations"), tabName = "information", icon = icon("info-circle", "fa-2x"))            
        )
    ),
    
    dashboardBody(
        shinyDashboardThemes(
            theme = "grey_light"
        ),
        
        tags$head( 
            tags$style(HTML(".main-sidebar { font-size: 20px; }")) #change the font size to 20
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
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/dag_2pops.pdf.png align="middle" height="auto" width="1275" margin="0 auto">'
        )
      }
    )
    
    output$models_picture <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/models.png align="middle" height="auto" width="1275" margin="0 auto">'
          )
        }
      )
    
    output$model_comparisons <-
      renderText({
        c(
          '<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/model_comparisons.png align="middle" height="auto" width="1275" margin="0 auto">'
        )
      }
      )
    
    
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
        a <-column(width = 10,
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
            a <- infoBox(title= NULL, value = h3("NON CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_upload == TRUE){
            a <- infoBox(title= NULL, value = h3("CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check filtering
    output$check_filtering_info <- renderUI({
        if(input$check_filtering == FALSE) {
            a <- infoBox(title= NULL, value = h3("NON CHECKED"), subtitle = NULL, icon = icon("bath"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_filtering == TRUE){
            a <- infoBox(title= NULL, value = h3("CHECKED"), subtitle = NULL, icon = icon("bath"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check populations
    output$check_populations_info <- renderUI({
        if(input$check_populations == FALSE) {
            a <- infoBox(title= NULL, value = h3("NON CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_populations == TRUE){
            a <- infoBox(title= NULL, value = h3("CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "green", fill = TRUE, width = 3)
        }
    })
    
    ## Check prior
    output$check_prior_info <- renderUI({
        if(input$check_prior == FALSE) {
            a <- infoBox(title= NULL, value = h3("NON CHECKED"), subtitle = NULL, icon = icon("dice"), color = "red", fill = TRUE, width = 3)
        } else if(input$check_prior == TRUE){
            a <- infoBox(title= NULL, value = h3("CHECKED"), subtitle = NULL, icon = icon("dice"), color = "green", fill = TRUE, width = 3)
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

