setwd("/Users/martinastoycheva/Desktop/Arabidopsis_iGA")
source('Francesco_code_modified.R')
library(Matrix)
library(ChemmineOB)
library(ChemmineR)
library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("METABOLiGA",
             h4("Metabolite iterative Group Analysis")),
  br(),
  br(),
  tabsetPanel(
    tabPanel("Overview",
             h3("Iterative Group Analysis (iGA)"),
             br(),
             h5("What is the idea behind iGA?"),
             div("iGA is a highly useful tool for biologist who want easy automatic functional analysis of"),
             div("'omics' experimental data. It is based on the concept that coordinated change of gene expression"), 
             div("or molecular concentration is physilogically relevent and performs functional set analysis of the data."),
             div("iGA was initially designed for microarray data analysis but is applicable to all omics data: genomics, proteomics,"),
             div("metabolomics. If molecules decrease and increase in the same time it is highly likely they are part of the same"),
             div("biological pathway or interrelated pathways. Now, there are two main questions that arise from this idea. 'What are"),
             div("the biological pathways involved?' and 'How likely is this?'"),
             div("iGA answers both of these questions by providing automatic group annotation that can be any"),
             div("available information about a molecule or a gene that the researcher has or can find, and by providing"),
             div("p value to show how significant the concerted change of all or some members of the group is."),
             br(),
             h5("How does iGA work?"),
             div("First, all of the genes, metabolites or proteins are arranged according to a metric chosen by the user."), 
             div("It can be t-statistic, f-statistic, fold change or anything considered relevant. iGA determines optimal"),
             div("thresholds for each group of genes, counting group members and checking through the list on each iteration."),
             div("The p-value is determined by the question: 'How likely is it to observe that many members of the class that"), 
             div("high up in the list by chance?' Third, the member of each class with the lowest p-value is found and this value"),
             div("is called the PC value. Last, the list is sorted by PC value (Breitling et al., 2004)."),
             br(),
             h3("How does METABOLiGA work?"),
             div("The web tool delivers automatic functional anotation for metabolomics experimental results using iGA function."),
             div("METABOLiGA is a biologists friendly user interface that provides easy to use iGA with additional functionalities  "),
             div("such as use of SMILES strings and PubChem IDs. The annotation is based on statistical calculation that detects "),
             div("similutataneous changes in concentration in either direction (decreasing or increasing). It provides a direct "), 
             div("answer with statistical significance - min.PC value i.e p value. It determines these functional groups from a"),
             div("list of metabolites that are enriched according to a metric which can be either their molecular weight calculated"), 
             div("from structure inferred from SMILES string, or to a metric provided by the user which can be any thing the user  "),
             div("considers relevant for instance differential expression. The metabolites can be a list of PubChem IDs, PubChem SMILES "),
             div("file as downloaded from the database, SMILES only or compound names. The app urges the user to input the correct"), 
             div("content with interactive conditional panels and multiple explanations of usage. Moreover, it produces easy to read "),
             div("output in the form of tables that are colourcode and can be searched and sorted by significance."),
             br(),
             br(),
             div("Breitling R., Amtmann A., Herzyk P. (2004). Iterative Group Analysis (iGA): A simple tool to enhance sensitivity and facilitate"), 
             div("interpretation of microarray experiments. BMC Bioinformatics, 5:34 "),
             br(),
             br(),
             br()
    ),
    tabPanel('File Upload',
             sidebarPanel(
               h2("Upload & Choose Dataset"),

               selectInput(
                 "met_label", "Choose a type of metabolite indicator input type",
                 c("Annotations User Input" = "annotations",
                   "PubChem IDs" = "PubChem_IDs",
                   "PubChem SMILES file" = "PubChemSMILES",
                   "SMILES" = "SMILES")),
               # Only show this panel if annotations file is chosen
               conditionalPanel(
                 condition = "input.met_label == 'annotations'",
                 # Input: Select a PubChem Smiles file
                 fileInput("file1", "Upload CSV annotations file",
                           multiple = FALSE) 
               ),
               # Only show this panel if PubChem Smiles file is chosen
               conditionalPanel(
                 condition = "input.met_label == 'PubChemSMILES'",
                 # Input: Select a PubChem Smiles file
                 fileInput("file2", "Upload TXT PubChem SMILES file",
                           multiple = FALSE) 
               ),
               # Only show this panel if SMILES file is chosen
               conditionalPanel(
                 condition = "input.met_label == 'SMILES'",
                 # Input: Select a PubChem Smiles file
               fileInput("file3", "Upload TXT SMILE strings file",
                         multiple = FALSE)
               ),
               # Only show this panel if the pub chem id file is chosen
               conditionalPanel(
                 condition = "input.met_label == 'PubChem_IDs'",
                 # Input: Select a PubChem IDs file
                 fileInput("file4", "Upload TXT PubChem IDs file",
                           multiple = FALSE) 
                ),
                
               selectInput(
                 "metric_type", "Choose a type of metric",
                 c("Molecular Weight" = "weight_metric",
                   "File Input Metric" = "metric_file")),
                conditionalPanel(
                  condition = "input.metric_type == 'metric_file'",
                  # Input: Select a Metric file
                  fileInput("file5", "Upload CSV Metric file",
                            multiple = FALSE),
                  # Input: Checkbox if file has header ----
                  checkboxInput("header", "Header", TRUE)
                 ),
               # Horizontal line ----
               tags$hr(),
               
               #Submit button to do t-test after data is uploaded
               actionButton(
                 inputId = "submit_iGA",
                 label = "Submit"),
               helpText("Note: The program will take a few minutes to run before displaying results in the results tab.")
             ),
             
             mainPanel(
               h2("Input Instructions"),
               br(),
               h3("Identifiers: "),
               div("The program requires a file that contains identifications for the molecules that will be analysed."),
               div("Four different files are accepted: "),
               withTags({
                 ol(
                   li("File (.csv) containing ready table: metabolite names/IDs of any database in row 1 and ready annotations for each"),
                   li("File (.txt) containing PubChem IDs i.e 6353"),
                   li("File (.txt) containing SMILES strings and PubChem IDs such as the output file downloaded from PubChem"),
                   li("File (.txt) containing only SMILES string")
                 )
               }),
                 
                 h3("Metric: "),
                 div("The program allows you to use the Molecular Weight calculated from the SMILES strings or to input a metric such as:"), 
                 div("fold change, f or t statistic. Note: metric file must be input if ready table option is chosen for identifiers! "),
                 div("The file has to be in .csv format and to have two columns: "),
                 withTags({
                   ol(
                     li("Column 1: Names or identifiers of the molecules."),
                     li("Column 2: The metric you would like to use.")
                   )
             }),
                div("The interface allows you to select: "),
               withTags({
                 ol(
                     li("Whether the first line in the metric file is a header or should be used.")
                 )
               }),
             div("  *Default is TRUE.")
             )
             
    ),
    tabPanel('Result Tables',
            br(),
            h3("iGA Summary Tables"),
            br(),
            h4("iGA Summary Table Decreasing"),
            DT::dataTableOutput("iga_summary"),
            br(),
            h4("iGA Summary Table Increasing"),
            DT::dataTableOutput("iga_summary2"),
            br(),
            tags$hr(),
            h3("Additional Information"),
            br(),
            h4("Total Group Number: "),
            verbatimTextOutput("groups.total"),
            br(),
            h4("iGA Results Information:"),
            h5("Deacreasing order iGA. PC values per group:"),
             verbatimTextOutput("groups"),
            h5("Increasing order iGA. PC values per group:"),
             verbatimTextOutput("groups2")
            
    )
    
    
    )
)

server <- function(input, output, session){
  
  #Increase max upload size from 5Mb (default) to 30Mb
  options(shiny.maxRequestSize=30*1024^2)
  
  #Data files secure upload
  data1 <- reactive({
    inFile1 <- input$file1
    if (is.null(inFile1)){return (NULL)}
    read.csv(input$file1$datapath)
  })
  
  #Data files secure upload
  data2 <- reactive({
    inFile1 <- input$file2
    if (is.null(inFile2)){return (NULL)}
    read.table(input$file2$datapath)
  })

  data3 <- reactive({
    inFile2 <- input$file3
    if (is.null(inFile3)){return (NULL)}
    read.table(input$file3$datapath)
  })

  data4 <- reactive({
    inFile3 <- input$file4
    if (is.null(inFile4)){return (NULL)}
    read.table(input$file4$datapath)
  })

  data5 <- reactive({
    inFile4 <- input$file5
    if (is.null(inFile5)){return (NULL)}
    read.csv(input$file5$datapath)
  })
  
  observeEvent(
    eventExpr = input$submit_iGA,
    handlerExpr = {
       
      #Load files for smiles, ids or annotations inside reactive event
      if (is.null(input$file1)){
        upload_annotations <- NULL  
      } else {
        upload_annotations <- read.csv(input$file1$datapath, header = F, row.names = 1)
      }  
      if (is.null(input$file2)){
        upload_pubchemsmiles <- NULL  
      } else {
        upload_pubchemsmiles <- read.table(input$file2$datapath)
      }
      if (is.null(input$file3)){
        upload_smiles <- NULL  
      } else {
        upload_smiles <- read.table(input$file3$datapath)
      }
      if (is.null(input$file4)){
        upload_pubchemids <- NULL  
      } else {
        upload_pubchemids <- read.table(input$file4$datapath)
      }
    # Conditional panel choice of file input and create iGA inputs
    if (input$met_label == 'annotations'){
      trp_matrix <- data.frame(upload_annotations)
      #Transpose to remove automatic column names
      trp_matrix <- t(trp_matrix)
      row.names(trp_matrix) <- NULL
      trp_matrix <- t(trp_matrix)
      #delete NA no information rows 
      trp_matrix <- trp_matrix[complete.cases(trp_matrix),]
      #Remove column names 
      remove_names <- trp_matrix
      row.names(remove_names) <- NULL
      remove_names <- remove_names[complete.cases(remove_names)]
      unique_annotations <- unique(remove_names)
      #Create empty matrix
      matrix_sample_membership <- Matrix(0,nrow(trp_matrix), length(unique_annotations))
      colnames(matrix_sample_membership) <- unique_annotations
      rownames(matrix_sample_membership) <- row.names(trp_matrix)
      #Now the prorgram checks whether unique_annotations matches anything in each row i.e each molecule in trp_matrix_KEGG
      for(i in 1:nrow(trp_matrix)){
        index <- which(unique_annotations%in%trp_matrix[[i]])
        matrix_sample_membership[i,index] <- rep(1, length(index))
      }
      names <- rownames(matrix_sample_membership)
    }
    else if (input$met_label == 'PubChemSMILES'){
      try_smiles <- upload_pubchemsmiles
      #Swap columns of PubChem SMILE file
      try_smiles <- try_smiles[c(2,1)]
      try_smiles <- as.matrix.data.frame(try_smiles)
      #write to file with correct format accepted by read.SMIset
      write.table(try_smiles, file = "try_smi.txt", sep = "\t", row.names = F, col.names = F,quote = F)
      try_smiles_smi <- read.SMIset("try_smi.txt")
      #Convert smiles data to sdf
      sdfset <- smiles2sdf(try_smiles_smi)
      names <- try_smiles[,2]
      #Make functional groups out of sdfset
      sample_fctgroup <- groups(sdfset, groups="fctgroup", type = "countMA")
      #Columns should not have 0 sum so are omitted
      my.mat <- replace(sample_fctgroup, sample_fctgroup > 1, 1)
      #Rename matrix 
      matrix_sample_membership <- Matrix(my.mat)
    }
     else if(input$met_label =='SMILES'){
       #SMILES only
       try_smiles <- upload_smiles
       try_smiles <- as.matrix.data.frame(try_smiles)
       #write to file with correct format accepted by read.SMIset
       write.table(try_smiles, file = "try_smi.txt", sep = "\t", row.names = F, col.names = F,quote = F)
       try_smiles_smi <- read.SMIset("try_smi.txt")
       sdfset <- smiles2sdf(try_smiles_smi)
       names <- NULL 
       #Make functional groups out of sdfset
       sample_fctgroup <- groups(sdfset, groups="fctgroup", type = "countMA")
       #Columns should not have 0 sum so are omitted
       my.mat <- replace(sample_fctgroup, sample_fctgroup > 1, 1)
       #Rename matrix 
       matrix_sample_membership <- Matrix(my.mat)
    } 
     else if(input$met_label =='PubChem_IDs'){ 
      upload_pubchemids <- unique(upload_pubchemids)
      try_ids <- as.vector(upload_pubchemids)
      try_ids <- as.numeric(unlist(try_ids))
      sdfset <- getIds(try_ids) 
      names <- sdfid(sdfset)
      #Make functional groups out of sdfset
      sample_fctgroup <- groups(sdfset, groups="fctgroup", type = "countMA")
      #Columns should not have 0 sum so are omitted
      my.mat <- replace(sample_fctgroup, sample_fctgroup > 1, 1)
      #Rename matrix 
      matrix_sample_membership <- Matrix(my.mat)
     }
    #Make metric 
    if (is.null(input$file5)){
      metric_upload <- NULL  
    }else{
      metric_upload <- read.csv(input$file5$datapath, header = input$header)
    }  

    if (input$metric_type == 'metric_file'){
    metric_use <- metric_upload
    metric_use <- metric_use[,2]
    }
    else if (input$metric_type == 'weight_metric'){
    metric_use <- as.matrix(MW(sdfset, addH = FALSE))
    metric_use <- metric_use[,1]
    }
    
    #Execute iGA: details on the input it accepts are available on the iGA github
    iGA_out <- iGA_acc(metric = metric_use, group.membership = matrix_sample_membership[, which(colSums(matrix_sample_membership)>0)],
                   groups = colnames(matrix_sample_membership[, which(colSums(matrix_sample_membership)>0)]), var.names = names, decreasing = TRUE)
    
    #Make a dataframe for the result table
    iGA_result_table <- data.frame(iGA_out$summary)
    #formating for scientific digits
    iGA_result_table <- format.data.frame(iGA_result_table, digits = 3)
    #output datatable for iGA summary
    if (is.null(iGA_result_table)){
      output$iga_summary <- NULL  
    }else{
      output$iga_summary <- DT::renderDataTable({
        datatable(iGA_result_table) %>%
          formatStyle('PC.adjusted',
                      backgroundColor = styleInterval(0.5, c('yellow', 'grey'))) %>%
          formatStyle('min.PC',
                      backgroundColor = styleInterval(0.5, c('yellow', 'grey')))
      }) 
    }
    
    #Execute iGA: details on the input it accepts are available on the iGA github
    iGA_out2 <- iGA_acc(metric = metric_use, group.membership = matrix_sample_membership[, which(colSums(matrix_sample_membership)>0)],
                       groups = colnames(matrix_sample_membership[, which(colSums(matrix_sample_membership)>0)]), var.names = names, decreasing = FALSE)
    
    #Make a dataframe for the result table
    iGA_result_table2 <- data.frame(iGA_out2$summary)  
    #formating for scientific digits
    iGA_result_table2 <- format.data.frame(iGA_result_table2, digits = 3)
    #output datatable for iGA summary
    if (is.null(iGA_result_table2)){
      output$iga_summary2 <- NULL  
    }else{
      output$iga_summary2 <- DT::renderDataTable({
        datatable(iGA_result_table2) %>%
          formatStyle('PC.adjusted',
                      backgroundColor = styleInterval(0.5, c('yellow', 'grey'))) %>%
          formatStyle('min.PC',
                      backgroundColor = styleInterval(0.5, c('yellow', 'grey')))
      }) 
    }    
    
    output$groups.total <- renderPrint({
      paste("Total number of groups: ", length(iGA_out$minPCs.pos))
    })
    
    output$groups <- renderPrint({
      paste(iGA_out$PC.list)
    })
    
    output$groups2 <- renderPrint({
      paste(iGA_out2$PC.list)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)