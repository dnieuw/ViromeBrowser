#GUI
importFilesOutput <- function(id){
  ns <- NS(id)
  tabBox(width = "100%",selected = "File Import",
    tabPanel(title = "File Import",
      fluidPage(
        column(width = 3,
          box(title = "File Import", solidHeader = T, status = "primary", width = NULL,
            box(width = 12,
              uiOutput(ns("metadata_file")),
              uiOutput(ns("annot_file")),
              uiOutput(ns("contig_file"))
            ),
            box(width = 12,
              selectInput(inputId = ns("load.data.func"),
                label = "Select annotation type:",
                # Currently there is only a function for importing BLAST files
                choices = c("BLAST format"="blast.data"),
                selected = "blast.data"
              ),
              actionButton(ns("load_data"), "Load Data"),
              br(),br(),
              prettySwitch(
                inputId = ns("load_example"),
                label = "Load example data", 
                status = "success",
                fill = TRUE
              )
            )
          )
        ),
        column(width = 9,
          importMetadataOutput(ns("metadata"))
        )
      )
    ),
    tabPanel(title = "Annotation Files",
      importAnnotationOutput(ns("annotation"))
    ),
    tabPanel(title = "Contig Files",
      importContigOutput(ns("contig"))
    )
  )
}

#callModule function
importFiles <- function(input, output, session){
  #Initialize namespace
  ns <- session$ns
  
  #Input annotation file GUI
  output$annot_file <- renderUI({
    fileInput(ns("annot"), 
              label = "Annotation files (.csv/.tsv)", 
              accept = c(".csv",".tsv"), 
              multiple = T)
  })
  
  #Input contig file GUI
  output$contig_file <- renderUI({
    fileInput(ns("contig"), 
              label = "Contig files (.fasta)", 
              accept = c(".fas",".fasta"), 
              multiple = T)
  })
  
  #Input metadata file GUI
  output$metadata_file <- renderUI({
    fileInput(ns("metadata"), 
              label = "Metadata file (.csv/.tsv)", 
              accept = c(".csv",".tsv"), 
              multiple = F)
  })
  
  annot <- callModule(importAnnotation,id = "annotation", reactive(input$annot), reactive(input$load.data.func), reactive(input$load_data), reactive(input$load_example))
  contig <- callModule(importContig, id = "contig", reactive(input$contig), reactive(input$load.data.func), reactive(input$load_data), reactive(input$load_example))
  metadata <- callModule(importMetadata, id = "metadata", reactive(input$metadata), reactive(input$load.data.func), reactive(input$load_data), reactive(input$load_example))
  
  return(list("annot" = annot$data, "defaults" = annot$defaults, "contig" = contig, "metadata" = metadata))
}

#GUI
importAnnotationOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Annotation files", solidHeader = T, status = "primary", width = NULL,
        DT::dataTableOutput(ns("example")),
        footer = "The first line of each imported and parsed file is shown"
    )
  )
}

#callModule function
importAnnotation <- function(input, output, session, annot_file, load.data.func, load.data, exampledata){
  
  #Reactive event for loading tsv data
  imported.data <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(annot_file())
    }
    
    return(get(load.data.func())())
  })
  
  #Custom data import function for slim output
  blast.data <- reactive({
    
    examplefiles <- function() {
      example_names <- c("example1.tsv","example2.tsv","example3.tsv","example4.tsv")
      example_file <- sapply(example_names, function(x) system.file("extdata", x, package = "viromeBrowser"))
      example_data <- data.frame("name"=example_names, "datapath"=example_file)
      return(example_data)
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading BLAST annotation files", value = 0)
    on.exit(progress$close())
    
    if(exampledata()==TRUE) {
      #Load the example file names if we want to load the example data
      file_data <- examplefiles()
    } else {
      file_data <- annot_file()
    }
    
    n <- nrow(file_data)
    
    data <- apply(file_data, 1, function(file){
      progress$inc(1/n, detail = paste0("Reading: ", file["name"]))
      result <- fread(file["datapath"])
      result$file.id <- str_remove(file["name"], "\\..*")
      return(result)
    })
    
    data <- rbindlist(data)

    names(data) <- c("contig.id","annotation","contig.length","hit.start","hit.end","length.homology","frac.homology","qframe", "identity","evalue","target","file.id")
    
    data[,contig.id := as.factor(contig.id)]
    data[,annotation := as.factor(annotation)]
    data[,contig.length := as.numeric(contig.length)]
    data[,hit.start := as.numeric(hit.start)]
    data[,hit.end := as.numeric(hit.end)]
    data[,length.homology := as.numeric(length.homology)]
    data[,frac.homology := as.numeric(frac.homology)]
    data[,qframe := as.factor(qframe)]
    data[,identity := as.numeric(identity)]
    data[,evalue := -log10(as.numeric(evalue))]
    data[is.infinite(evalue) & sign(evalue) == 1, evalue := 999]
    data[is.infinite(evalue) & sign(evalue) == 1, evalue := 0]
    data[,"evalue(-log10)" := round(evalue, 2)]
    data[,evalue:=NULL]
    data[,target := as.factor(target)]
    data[,file.id := as.factor(file.id)]
    
    data <- data[,c("file.id","contig.id","annotation","contig.length","hit.start","hit.end","length.homology","frac.homology","qframe","identity","evalue(-log10)","target")]

    return(data)
  })
  
  #Table output of first line of each imported file
  output$example <- DT::renderDataTable({
    req(imported.data())
    #Group data by file.id and show first entry per file
    DT::datatable(imported.data()[, .SD[1, ], by = .(file.id)],
              rownames = FALSE,
              options = list(
                dom = "tip",
                pageLength = 15,
                scrollX = TRUE,
                sScrollY = "650px"
              )
    )
  })
  
  default_threshold_values <- list("file.id"=NA, "contig.id"=NA, "annotation"=NA,
                         "contig.length"=c(500,Inf),"hit.start"=NA,"hit.end"=NA,
                         "length.homology"=c(500,Inf),"frac.homology"=NA,"identity"=c(90,Inf),
                         "evalue(-log10)"=NA,"target"=NA)
  
  #Return the imported data
  return(list("data" = imported.data, "defaults" = default_threshold_values))
}

#GUI
importContigOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Contig files", solidHeader = T, status = "primary", width = NULL,
        DT::dataTableOutput(ns("example")),
        footer = "The first contig of each imported and parsed file is shown"
    )
  )
}

#CallModule function
importContig <- function(input, output, session, contig_file, load.data.func, load.data, exampledata){
  
  #Reactive event for loading contig data
  imported.fasta <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(contig_file())
    }

    return(get(load.data.func())())
  })
  
  #Custom data import function for blast output
  blast.data <- reactive({
    
    examplefiles <- function() {
      example_names <- c("example1.fasta","example2.fasta","example3.fasta","example4.fasta")
      example_file <- sapply(example_names, function(x) system.file("extdata", x, package = "viromeBrowser"))
      example_data <- data.frame("name"=example_names, "datapath"=example_file)
      return(example_data)
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading contig files", value = 0)
    on.exit(progress$close())
    
    if(exampledata()==TRUE) {
      #Load the example file names if we want to load the example data
      file_data <- examplefiles()
    } else {
      file_data <- contig_file()
    }
    
    n <- nrow(file_data)
    
    data <- apply(file_data, 1, function(file){
      #Increase progress
      progress$inc(1/n, detail = paste0(file["name"],"-indexing"))
      
      #Open and index fasta file
      fa <- open(FaFile(file["datapath"]))
      idx <- scanFaIndex(fa)
      
      return(list("fasta"=fa,"idx"=idx))
    })
    names(data) <- str_remove(file_data$name, "\\..*")

    return(data)
  })
  
  output$example <- DT::renderDataTable({
    data <-	req(imported.fasta())
    printdata <- lapply(seq_along(data), function(n){
      contig.example <- scanFa(data[[n]]$fasta, param=data[[n]]$idx[1])
      file.id <- names(data)[n]
      contig.id <- names(contig.example)
      contig.seq <- paste(as.character(subseq(contig.example[[1]], end = 15)),as.character(subseq(contig.example[[1]], start = -15)),sep="...")
      data.table("file.id"=file.id, "contig.id"=contig.id, "contig.seq"=contig.seq)
    })
    printdata <- rbindlist(printdata)
    DT::datatable(printdata,
              rownames = FALSE,
              options = list(
                dom = "tip",
                pageLength = 15,
                scrollX = TRUE,
                sScrollY = "650px"
              )
    )
  })
  
  return(imported.fasta)
}

#GUI
importMetadataOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Metadata", solidHeader = T, status = "primary", width = NULL,
        DT::dataTableOutput(ns("example"))
    )
  )
}

#CallModule function
importMetadata <- function(input, output, session, metadata_file, load.data.func, load.data, exampledata){
  
  #Reactive event for loading metadata data
  imported.metadata <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(metadata_file())
    }
    
    return(get(load.data.func())())
  })
  
  #Custom data import function for blast output
  blast.data <- reactive({
    
    examplefiles <- function() {
      example_name <- "example_metadata.txt"
      example_file <- system.file("extdata", "example_metadata.tsv", package = "viromeBrowser")
      example_data <- data.frame("name"=example_name, "datapath"=example_file)
      return(example_data)
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading metadata file", value = 0)
    on.exit(progress$close())
    
    if(exampledata()==TRUE) {
      #Load the example file names if we want to load the example data
      file_data <- examplefiles()
    } else {
      file_data <- metadata_file()
    }
    
    data <- fread(as.character(file_data[,"datapath"]))
    
    #First column should be the file.id key value
    names(data)[1] <- "file.id"

    #Remove all NA columns
    data <- data[,data[ ,sapply(.SD, function(x) !all(is.na(x)))], with=F]
    
    #Remove all columns with only one value
    data <- data[,data[ ,sapply(.SD, function(x) length(unique(x)) > 1)], with=F]

    #Look for dates and transform them to dates
    can_convert <- names(data)[which(!is.na(as.Date(unlist(data[1,]),"%m/%d/%Y")))]
    replacement_dates <- lapply(can_convert, function(x) data[,as.Date(get(x), "%m/%d/%Y")])
    data[,(can_convert) := replacement_dates]
    
    #Turn all character columns into factors
    character_columns <- names(data)[which(data[ ,sapply(.SD, class)]=="character")]
    replacement_factors <- data[,character_columns, with=F][,lapply(.SD, as.factor)]
    data[,(character_columns) := replacement_factors]
    
    return(data)
  })
  
  output$example <- DT::renderDataTable({
    tabledata <-	req(imported.metadata())
    
    DT::datatable(tabledata,
                  rownames = FALSE,
                  options = list(
                    dom = "tip",
                    pageLength = 15,
                    scrollX = TRUE,
                    sScrollY = "650px"
                  )
    )
  })
  
  return(imported.metadata)
}

annotationTableOutput <- function(id) {
  ns <- NS(id)
  fluidPage(
    box(width = 12,
        DT::dataTableOutput(ns("table")),
        actionButton(ns("clearbutton"), "Clear Selection")
    )
  )
}

annotationTable <- function(input, output, session, datasubset) {
  ns <- session$ns
  
  #Generate a table with contig info
  output$table <- DT::renderDataTable({
    DT::datatable(datasubset(),
                  selection = list(
                    mode = 'multiple',
                    selected = NULL,
                    target = 'row'
                  ),
                  options = list(
                    ordering=F,
                    dom = "tip",
                    pageLength = 100,
                    scrollX = TRUE,
                    scrollY = 100,
                    sScrollY = "500px"
                  )
    )
  })
  
  proxy = DT::dataTableProxy('table')
  
  observeEvent(input$clearbutton, {
    proxy %>% DT::selectRows(NULL)
  })
  
  get.selected <- reactive({
    return(input$table_rows_selected)
  })
  
  return(list("table" = datasubset, "selected" = get.selected))
}

