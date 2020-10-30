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
  
  annot <- callModule(importAnnotation,id = "annotation", reactive(input$annot), reactive(input$load_data), reactive(input$load_example))
  contig <- callModule(importContig, id = "contig", reactive(input$contig), reactive(input$load_data), reactive(input$load_example))
  metadata <- callModule(importMetadata, id = "metadata", reactive(input$metadata), reactive(input$load_data), reactive(input$load_example))
  
  return(list("annot" = annot$data, "defaults" = annot$defaults, "contig" = contig, "metadata" = metadata))
}

#GUI
importAnnotationOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    column(4,
      box(title = "Annotation configuration", solidHeader = T, status = "primary", width = NULL,
        rHandsontableOutput(ns("import_config_table"), width="100%"),
        footer = "Please select which columns correspont to which annotations in the table"
      )
    ),
    column(8,
      box(title = "Annotation files", solidHeader = T, status = "primary", width = NULL,
        dataTableOutput(ns("example")),
        footer = "The first line of each imported and parsed file is shown"
      )
    )
  )
}

#callModule function
importAnnotation <- function(input, output, session, annot_file, load.data, exampledata){
  
  #Reactive event for loading tsv data
  imported.data <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(annot_file())
    } 
    
    examplefiles <- function() {
      example_names <- c("example1.tsv","example2.tsv","example3.tsv","example4.tsv")
      example_file <- sapply(example_names, function(x) system.file("extdata", x, package = "viromeBrowser"))
      example_data <- data.frame("name"=example_names, "datapath"=example_file)
      return(example_data)
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading annotation files", value = 0)
    on.exit(progress$close())
    
    if (exampledata() == FALSE){
      file_data <- annot_file()
    } else {
      file_data <- examplefiles()
    }

    n <- nrow(file_data)
    
    data <- apply(file_data, 1, function(file){
      progress$inc(1/n, detail = paste0("Reading: ", file["name"]))
      result <- fread(file["datapath"])
      result$file.id <- str_remove(file["name"], "\\..*")
      return(result)
    })
    
    data <- rbindlist(data)
    
    return(data)
  })
  
  #Generate a table to select which columns represent what data to import
  output$import_config_table <- renderRHandsontable({
    req(imported.data())
    
    #Create a dataframe with all annotation data that can be imported and the column names from the imported data
    #Users can select which column in the input contains which annoation data
    input_options <- factor("NotAvailable",levels=c("NotAvailable", colnames(imported.data())))
    selection_df <- data.frame(column=rep(input_options, 17), 
                               row.names = c("contig.id",
                                    "superkingdom",
                                    "kingdom",
                                    "class",
                                    "order",
                                    "family",
                                    "genus",
                                    "species",
                                    "taxonomic_name",
                                    "contig.length",
                                    "hit.start",
                                    "hit.end",
                                    "length.homology",
                                    "frac.homology",
                                    "identity",
                                    "evalue",
                                    "target"))
    
    rhandsontable(selection_df,
                  colHeaders = FALSE,
                  rowHeaderWidth = 150,
                  stretchH = "all"
                  )
  })
  
  output$example <- renderDataTable({
    req(configured_data())
    
    #Group data by file.id and show first entry per file
    datatable(configured_data()[,.SD[1], by=file.id],
                  rownames = FALSE,
                  options = list(
                    dom = "tip",
                    pageLength = 15,
                    scrollX = TRUE,
                    sScrollY = "650px"
                  )
    )
    
  })
  
  configured_data <- reactive({
    
    if (exampledata() == FALSE){
      req(input$import_config_table)
      config <- as.data.table(hot_to_r(input$import_config_table), keep.rownames = T)
      colnames(config) <- c("rn","column")
      
      #Dont show anything if all columns are excluded
      if (config[,all(column=="NotAvailable")]) {
        return(NULL)
      }
      
      input_columns <- config[column!="NotAvailable",as.character(column)]
      annotation_config <- config[column!="NotAvailable",rn]
      
      data <- imported.data()
      data <- data[,c("file.id",input_columns), with=F]
      
      colnames(data) <- c("file.id", annotation_config)
    } else {
      data <- imported.data()
    }

    #Reformat the annotation columns if present
    sapply(colnames(data), switch,
      contig.id = {
        data[,contig.id := as.factor(contig.id)]
      },
      superkingdom = {
        data[,superkingdom := as.factor(superkingdom)]
      },
      kingdom = {
        data[,kingdom := as.factor(kingdom)]
      },
      class = {
        data[,class := as.factor(class)]
      },
      order = {
        data[,order := as.factor(order)]
      },
      family = {
        data[,family := as.factor(family)]
      },
      genus = {
        data[,genus := as.factor(genus)]
      },
      species = {
        data[,species := as.factor(species)]
      },
      taxonomic_name = {
        data[,taxonomic_name := as.factor(taxonomic_name)]
      },
      contig.length = {
        data[,contig.length := as.numeric(contig.length)]
      },
      hit.start = {
        data[,hit.start := as.numeric(hit.start)]
      },
      hit.end = {
        data[,hit.end := as.numeric(hit.end)]
      },
      length.homology = {
        data[,length.homology := as.numeric(length.homology)]
      },
      frac.homology = {
        data[,frac.homology := as.numeric(frac.homology)]
      },
      identity = {
        data[,identity := as.numeric(identity)]
      },
      evalue = {
        data[,evalue := -log10(as.numeric(evalue))]
        data[is.infinite(evalue) & sign(evalue) == 1, evalue := 999]
        data[is.infinite(evalue) & sign(evalue) == 1, evalue := 0]
        data[,"evalue(-log10)" := round(evalue, 2)]
        data[,evalue:=NULL]
      },
      target = {
        data[,target := as.factor(target)]
      },
      file.id = {
        data[,file.id := as.factor(file.id)]
      }
    )
    
    return(data)
    
  })
  
  default_threshold_values <- reactive({

    req(configured_data())
    
    thresholds <- NULL
    
    #Set the default thresholds for annotation columns if present
    sapply(colnames(data), switch,
           contig.id = {
             thresholds["file.id"]=NA
           },
           contig.length = {
             thresholds["contig.length"]=c(500,Inf)
           },
           length.homology = {
             thresholds["length.homology"]=c(500,Inf)
           },
           frac.homology = {
             thresholds["frac.homology"]=NA
           },
           identity = {
             thresholds["identity"]=c(90,Inf)
           },
           "evalue(-log10)" = {
             thresholds["evalue(-log10)"]=NA
           },
           target = {
             thresholds["target"]=NA
           },
           file.id = {
             thresholds["file.id"]=NA
           })
    
    return(thresholds)
    
  })
  
  #Return the imported data
  return(list("data" = configured_data, "defaults" = default_threshold_values))
}

#GUI
importContigOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Contig files", solidHeader = T, status = "primary", width = NULL,
        dataTableOutput(ns("example")),
        footer = "The first contig of each imported and parsed file is shown"
    )
  )
}

#CallModule function
importContig <- function(input, output, session, contig_file, load.data, exampledata){
  
  #Reactive event for loading contig data
  imported.fasta <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(contig_file())
    }
    
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
  
  output$example <- renderDataTable({
    data <-	req(imported.fasta())
    printdata <- lapply(seq_along(data), function(n){
      contig.example <- scanFa(data[[n]]$fasta, param=data[[n]]$idx[1])
      file.id <- names(data)[n]
      contig.id <- names(contig.example)
      contig.seq <- paste(as.character(subseq(contig.example[[1]], end = 15)),as.character(subseq(contig.example[[1]], start = -15)),sep="...")
      data.table("file.id"=file.id, "contig.id"=contig.id, "contig.seq"=contig.seq)
    })
    printdata <- rbindlist(printdata)
    datatable(printdata,
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
        dataTableOutput(ns("example"))
    )
  )
}

#CallModule function
importMetadata <- function(input, output, session, metadata_file, load.data, exampledata){
  
  #Reactive event for loading metadata data
  imported.metadata <- eventReactive(load.data(),{
    if (exampledata() == FALSE){
      req(metadata_file())
    }
    
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
  
  output$example <- renderDataTable({
    tabledata <-	req(imported.metadata())
    
    datatable(tabledata,
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
        dataTableOutput(ns("table")),
        actionButton(ns("clearbutton"), "Clear Selection")
    )
  )
}

annotationTable <- function(input, output, session, datasubset) {
  ns <- session$ns
  
  #Generate a table with contig info
  output$table <- renderDataTable({
    datatable(datasubset(),
                  selection = list(
                    mode = 'multiple',
                    selected = NULL,
                    target = 'row'
                  ),
                  options = list(
                    ordering=F,
                    dom = "tip",
                    pageLength = 10,
                    scrollX = TRUE
                    )
    )
  })
  
  proxy = dataTableProxy('table')
  
  observeEvent(input$clearbutton, {
    proxy %>% selectRows(NULL)
  })
  
  get.selected <- reactive({
    return(input$table_rows_selected)
  })
  
  return(list("table" = datasubset, "selected" = get.selected))
}

