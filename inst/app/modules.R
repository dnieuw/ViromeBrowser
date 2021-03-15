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
              uiOutput(ns("contig_file")),
              uiOutput(ns("bam_file"))
            ),
            box(width = 12,
              actionButton(ns("load_data"), "Load Data")
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
    ),
    tabPanel(title = "BAM Files",
      importBAMOutput(ns("bam"))
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
  
  #Input bam file GUI
  output$bam_file <- renderUI({
    fileInput(ns("bam"), 
              label = "BAM files (.bam)", 
              accept = c(".bam"), 
              multiple = T)
  })
  
  #Input metadata file GUI
  output$metadata_file <- renderUI({
    fileInput(ns("metadata"), 
              label = "Metadata file (.csv/.tsv)", 
              accept = c(".csv",".tsv"), 
              multiple = F)
  })
  
  annot <- callModule(importAnnotation,id = "annotation", reactive(input$annot), reactive(input$load_data))
  contig <- callModule(importContig, id = "contig", reactive(input$contig), reactive(input$load_data))
  readcounts <- callModule(importBAM, id = "bam", reactive(input$bam), reactive(input$load_data))
  metadata <- callModule(importMetadata, id = "metadata", reactive(input$metadata), reactive(input$load_data))
  
  return(list("annot" = annot$data, "defaults" = annot$defaults, "contig" = contig, "readcounts" = readcounts, "metadata" = metadata))
}

#GUI
importAnnotationOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Annotation files", solidHeader = T, status = "primary", width = NULL,
      uiOutput(ns("select_file")),
      DT::dataTableOutput(ns("example"))
    )
  )
}

#callModule function
importAnnotation <- function(input, output, session, annot_file, load.data){
  ns <- session$ns
  
  #Reactive event for loading tsv data
  imported.data <- eventReactive(load.data(),{
    
    req(annot_file())
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading annotation files", value = 0)
    on.exit(progress$close())
    
    file_data <- annot_file()

    n <- nrow(file_data)
    
    data <- apply(file_data, 1, function(file){
      progress$inc(1/n, detail = paste0("Reading: ", file["name"]))
      result <- fread(file["datapath"])
      names(result) <- c("contig.id","target","identity",
                       "length.homology","mismatch","gapopen",
                       "qstart","qend","hit.start","hit.end",
                       "evalue","bitscore","annotation")
      result$file.id <- str_remove(file["name"], "\\..*")
      return(result)
    })
    
    data <- rbindlist(data)
    
    return(data)
  })
  
  taxdump_lineage <- reactive({

    cache_dir <- tempdir()
    cache_file <- paste(cache_dir,"taxdump_tmp",sep='/')
    
    if (!file.exists(cache_file)) {
      tryCatch({
        progress <- shiny::Progress$new()
        progress$set(message = "Downloading taxdump", value = 0)
        on.exit(progress$close())
        download.file("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip",cache_file)
        progress$inc(1/2, detail = "Unzipping..")
        unzip(cache_file, exdir=cache_dir)
        progress$inc(1/2)
      },error=function(cond) {
        message("Cannot download new_taxdump")
        message(cond)
      }
      )
    }
    
    taxdata <- fread(file = paste(cache_dir,"rankedlineage.dmp",sep='/'),sep="\t")[,c(1,3,5,7,9,11,13,15,17,19)]
    colnames(taxdata) <- c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom")
    
    taxdata[,tax_id:=as.integer(tax_id)]
    setkey(taxdata,"tax_id")
    return(taxdata)
  })
  
  resolvemultipletaxon <- function(multitaxon, taxdump_lineage) {
    
    #Get entries with multiple taxon information and paste together all taxons per contig
    lineage_data <- multitaxon[,.(combined=paste(annotation,collapse = ";")),.(contig.id,file.id)]
    
    #Split on ";" to get 1 line per taxon id
    lineage_data <- lineage_data[,str_split(combined,";"),.(contig.id,file.id)]
    colnames(lineage_data) <- c("contig.id","file.id","tax_id")
    lineage_data[,tax_id:=as.integer(tax_id)]
  
    #Merge with ranked lineage to get complete lineage
    lineage_data <- merge.data.table(lineage_data,taxdump_lineage(), by = "tax_id")
    
    #Remove taxid which interferes with further process
    lineage_data[,tax_id:=NULL]
    
    find_uniq <- function(x) {
      if (all(x=="")) return("")
      x_uniq <- unique(x[x!=""])
      ifelse(length(x_uniq)==1,x_uniq,"-")
    }
    
    #Determine LCA per contig
    contig_lca <- lineage_data[, lapply(.SD,find_uniq), .(contig.id, file.id)]
    
    #Remove lineages with nonconforming lower ancestor
    contig_lca[superkingdom=='-',superkingdom:="root"]
    contig_lca[superkingdom=='-',c("tax_name","species","genus","family","order","class","phylum","kingdom"):=""]
    contig_lca[kingdom=='-',c("tax_name","species","genus","family","order","class","phylum","kingdom"):=""]
    contig_lca[phylum=='-',c("tax_name","species","genus","family","order","class","phylum"):=""]
    contig_lca[class=='-',c("tax_name","species","genus","family","order","class"):=""]
    contig_lca[order=='-',c("tax_name","species","genus","family","order"):=""]
    contig_lca[family=='-',c("tax_name","species","genus","family"):=""]
    contig_lca[genus=='-',c("tax_name","species","genus"):=""]
    contig_lca[species=='-',c("tax_name","species"):=""]
    contig_lca[tax_name=='-',c("tax_name"):=""]
    
    result <- merge.data.table(multitaxon, contig_lca, by=c("file.id","contig.id"))
    
    return(result)
  }
  
  configured.data <- reactive({
    req(imported.data())
    
    data <- imported.data()
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Adding taxonomic lineage to annotation", value = 0)
    progress$inc(0.5)
    on.exit(progress$close())
      
    req(taxdump_lineage())
    
    data[,result.id:=seq(nrow(data))]
    
    #Get entries with mutiple taxIDs
    multitaxon <- data[grep(";",annotation)]
    #Resolve results with multiple taxon ids using LCA
    multitaxon <- resolvemultipletaxon(multitaxon, taxdump_lineage)
    
    #Get entries with single taxon
    singletaxon <- data[grep(";",annotation,invert = T)]
    singletaxon[,annotation:=as.integer(annotation)]
    #Merge with ranked lineage to get complete lineage
    singletaxon <- merge.data.table(singletaxon,taxdump_lineage(), by.x = "annotation", by.y = "tax_id")
    
    data <- rbind.data.frame(singletaxon, multitaxon)
    progress$inc(0.5)
    
    #Remove taxID which is replaced by lineage information
    data[,annotation := NULL]
    
    data[,contig.id := as.factor(contig.id)]
    data[,hit.start := as.numeric(hit.start)]
    data[,hit.end := as.numeric(hit.end)]
    data[,length.homology := as.numeric(length.homology)]
    data[,identity := as.numeric(identity)]
    data[,bitscore := as.numeric(bitscore)]
    data[,evalue := -log10(as.numeric(evalue))]
    data[is.infinite(evalue) & sign(evalue) == 1, evalue := 999]
    data[is.infinite(evalue) & sign(evalue) == 1, evalue := 0]
    data[,"evalue(-log10)" := round(evalue, 2)]
    data[,evalue:=NULL]
    data[,target := as.factor(target)]
    data[,file.id := as.factor(file.id)]
    
    data <- data[,c("file.id","contig.id","target","tax_name",
                    "species","genus","family","order",
                    "class","phylum","kingdom","superkingdom",
                    "hit.start","hit.end","length.homology",
                    "identity","evalue(-log10)", "bitscore")]
    
    return(data)
  })
  
  output$select_file <- renderUI({
    req(imported.data())
    filenames <- unique(imported.data()$file.id)
    selectInput(ns("select_file"),label = "File", choices = filenames, width = "200px", multiple = F)
  })
  
  #Create example to show in the interface
  output$example <- DT::renderDataTable({
    req(configured.data(), input$select_file)
    #Group data by file.id and show first entry per file
    DT::datatable(configured.data()[file.id==input$select_file],
                  rownames = FALSE,
                  options = list(
                    dom = "tip",
                    pageLength = 15,
                    scrollX = TRUE,
                    sScrollY = "650px"
                  )
    )
  })
  
  #Set the default quality threshold values
  default_threshold_values <- list("length.homology"=list(search='500 ... Inf'),"contig.length"=list(search='500 ... Inf'),"identity"=list(search='90 ... Inf'))
  
  #Return the imported data
  return(list("data" = configured.data, "defaults" = default_threshold_values))
}

#GUI
importContigOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "Contig files", solidHeader = T, status = "primary", width = NULL,
        uiOutput(ns("select_file")),
        DT::dataTableOutput(ns("example"))
    )
  )
}

#CallModule function
importContig <- function(input, output, session, contig_file, load.data){
  ns <- session$ns
  #Reactive event for loading contig data
  imported.fasta <- eventReactive(load.data(),{
   
    req(contig_file())

    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading contig files", value = 0)
    on.exit(progress$close())
    
    file_data <- contig_file()
    
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
  
  output$select_file <- renderUI({
    req(imported.fasta())
    
    filenames <- names(imported.fasta())
    selectInput(ns("select_file"),label = "File", choices = filenames, width = "200px", multiple = F)
  })
  
  output$example <- DT::renderDataTable({
    req(input$select_file)
    data <-	req(imported.fasta())
    
    contig.example <- scanFa(data[[input$select_file]]$fasta, param=data[[input$select_file]]$idx)
    file.id <- input$select_file
    contig.id <- names(contig.example)
    contig.seq <- paste(as.character(subseq(contig.example[[1]], end = 15)),as.character(subseq(contig.example[[1]], start = -15)),sep="...")
    printdata <- data.table("file.id"=file.id, "contig.id"=contig.id, "contig.seq"=contig.seq)

    DT::datatable(printdata[file.id==input$select_file],
              rownames = FALSE,
              options = list(
                dom = "tip",
                pageLength = 15,
                scrollX = TRUE,
                sScrollY = "650px",
                initComplete = JS(
                  "function(settings, json) {",
                  "$('td').css({'border': '1px blue'});",
                  "$('th').css({'border': '1px red'});",
                  "}")
              )
    )
  })
  
  return(imported.fasta)
}

#GUI
importBAMOutput <- function(id){
  ns <- NS(id)
  fluidPage(
    box(title = "BAM files", solidHeader = T, status = "primary", width = NULL,
        uiOutput(ns("select_file")),
        DT::dataTableOutput(ns("example"))
    )
  )
}

#CallModule function
importBAM <- function(input, output, session, BAM_file, load.data){
  ns <- session$ns
  #Reactive event for loading contig data
  imported.bam <- eventReactive(load.data(),{
    
    req(BAM_file())
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading BAM files", value = 0)
    on.exit(progress$close())
    
    file_data <- BAM_file()
    
    #2 progress steps per file
    n <- nrow(file_data)*2
    
    data <- apply(file_data, 1, function(file){
      #Increase progress
      progress$inc(1/n, detail = paste0(file["name"],"-indexing"))
      
      #Index bamfile
      bam_index <- indexBam(file["datapath"])
      
      #Increase progress
      progress$inc(1/n, detail = paste0(file["name"],"-counting"))
      
      #Open bamfile
      bam_file <- BamFile(file["datapath"])
        
      #Load contig information (name,length)
      contig_info <- seqinfo(bam_file)
      contig_info <- data.table(seqnames=seqnames(contig_info),seqlengths=seqlengths(contig_info))
      
      #Create a genomic range for each complete contig
      gr <- GRanges(contig_info$seqnames, IRanges(start=1,end=contig_info$seqlengths))
      
      #Generate the parameters to scan the bamfile with
      params <- ScanBamParam(which = gr, what = scanBamWhat())
      
      #Extract read counts per contig for the complete contig (and contig size)
      counts <- data.table(countBam(bam_file, param = params)[,c(1,6)])
      
      names(counts) <- c("contig.id","readcount")
      
      counts[,totalreadcount:=sum(readcount)]
      
      counts[,file.id:=str_remove(file["name"], "\\..*")]
      
      return(counts)
    })
    
    data <- rbindlist(data)
    
    return(data)
  })
  
  output$select_file <- renderUI({
    req(imported.bam())
    filenames <- unique(imported.bam()$file.id)
    selectInput(ns("select_file"),label = "File", choices = filenames, width = "200px", multiple = F)
  })
  
  output$example <- DT::renderDataTable({
    req(input$select_file)
    data <-	req(imported.bam())
    
    
    DT::datatable(imported.bam()[file.id==input$select_file],
              rownames = FALSE,
              options = list(
                dom = "tip",
                pageLength = 15,
                scrollX = TRUE,
                sScrollY = "650px"
              )
    )
  })
  
  return(imported.bam)
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
importMetadata <- function(input, output, session, metadata_file, load.data){
  
  #Reactive event for loading metadata data
  imported.metadata <- eventReactive(load.data(),{
    
    req(metadata_file())
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Loading metadata file", value = 0)
    on.exit(progress$close())
    
    file_data <- metadata_file()

    data <- fread(as.character(file_data[,"datapath"]))
    
    #First column should be the file.id key value
    names(data)[1] <- "file.id"

    #Remove all NA columns
    data <- data[,data[ ,sapply(.SD, function(x) !all(is.na(x)))], with=F]
    
    #Remove all columns with only one value
    data <- data[,data[ ,sapply(.SD, function(x) length(unique(x)) > 1)], with=F]

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
                    ordering=T,
                    dom = "tipf",
                    pageLength = 15,
                    scrollX = TRUE,
                    scrollY = TRUE
                    )
    )
  })
  
  get.selected <- reactive({
    return(input$table_rows_selected)
  })
  
  return(list("table" = datasubset, "selected" = get.selected))
}

