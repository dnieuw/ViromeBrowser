#GUI
importFilesOutput <- function(id){
  ns <- NS(id)
  tabBox(width = "100%",selected = "Analysis files",
    tabPanel(title = "Analysis files",
      fluidPage(
        column(width = 3,
          box(title = "Input files", solidHeader = T, status = "primary", width = NULL,
            box(title = "Annotation file", width = 12,
              uiOutput(ns("annot_file"))
            ),
            box(title = "Contig file", width = 12,
              uiOutput(ns("contig_file"))
            ),
            box(width = 12,
              selectInput(inputId = ns("load.data.func"),
                label = "Select data type:",
                # Currently there is only functionallity for importing BLAST files
                choices = c("BLAST data"="blast.data", "Other"="other.data"),
                selected = "blast.data"
              ),
              actionButton(ns("load_data"),"Load the data"),
              prettySwitch(
                inputId = ns("load_example"),
                label = "Example data", 
                status = "success",
                fill = TRUE
              )
            )
          )
        ),
        column(width = 9,
          box(title = "Annotation file format", solidHeader = T, status = "primary", width = NULL,
            uiOutput(ns("format_info"))
          )
        )
      )
    ),
    tabPanel(title = "Annotation Input",
      importAnnotationOutput(ns("annotation"))
    ),
    tabPanel(title = "Contig Input",
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
              label = "Contig annotation file", 
              accept = c(".csv",".tsv"), 
              multiple = T)
  })
  
  #Input contig file GUI
  output$contig_file <- renderUI({
    fileInput(ns("contig"), 
              label = "Contig sequence file", 
              accept = c(".fas",".fasta"), 
              multiple = T)
  })
  
  output$format_info <- renderUI({
    req(input$load.data.func)
    switch(input$load.data.func, 
      "blast.data" = img(src='Annotation_format.png', align = "middle", width="100%")
    )
  })
  
  annot <- callModule(importAnnotation,id = "annotation", reactive(input$annot), reactive(input$load.data.func), reactive(input$load_data), reactive(input$load_example))
  contig <- callModule(importContig, id = "contig", reactive(input$contig), reactive(input$load.data.func), reactive(input$load_data), reactive(input$load_example))
  
  return(list("annot" = annot, "contig" = contig))
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
      example_names <- c("example1.tsv","example2.tsv","example3.tsv","example4.tsv","example5.tsv")
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
      progress$inc(1/n, detail = file["name"])
      result <- fread(file["datapath"])
      result$file.id <- str_remove(file["name"], "\\..*")
      return(result)
    })
    
    data <- rbindlist(data)

    names(data) <- c("contig.id","annotation","contig.length","hit.start","hit.end","length.homology","frac.homology","qframe", "identity","evalue","target","file.id")
    class(data$contig.length) <- "numeric"
    class(data$hit.start) <- "numeric"
    class(data$hit.end) <- "numeric"
    class(data$length.homology) <- "numeric"
    class(data$frac.homology) <- "numeric"
    class(data$identity) <- "numeric"
    class(data$evalue) <- "numeric"

    data <- data[,c("file.id","contig.id","annotation","contig.length","hit.start","hit.end","length.homology","frac.homology","identity","evalue","target")]
    
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
  
  #Return the imported data
  return(imported.data)
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
      example_names <- c("example1.fasta","example2.fasta","example3.fasta","example4.fasta","example5.fasta")
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