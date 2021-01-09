output$metadata_filter_selection <- renderUI({
  pickerInput(
    inputId = "metadata_filter_selection",
    label = "Filter by:", 
    choices = names(metadata()),
    options =  list(
      size = 5,
      "live-search" = TRUE,
      "max-options" = 9,
      "max-options-text" = "Maximum filters reached"
    ), 
    multiple = TRUE
  )
})

output$metadata_stratify_selection <- renderUI({
  pickerInput(
    inputId = "metadata_stratify_selection",
    label = "Stratify by:", 
    choices = names(metadata()),
    options =  list(
      size = 5,
      "live-search" = TRUE
    ), 
    multiple = FALSE
  )
})

output$annotation_count_selection <- renderUI({

  choices <- c("Contig count"="countcontig",
              "Absolute readcount" = "abs_readcount",
              "Total scaled readcount" = "scaled_readcount")

  pickerInput(
    inputId = "annotation_count_selection",
    label = "Fill with:",
    choices = choices,
    options =  list(
      size = 5
    ),
    multiple = FALSE
  )
})

output$annotation_stratify_selection <- renderUI({

  annotation_levels <- c("superkingdom",
                         "kingdom",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "species",
                         "tax_name")

  pickerInput(
    inputId = "annotation_stratify_selection",
    label = "Annotate with:",
    choices = annotation_levels,
    selected = "family",
    options =  list(
      size = 5
    ),
    multiple = FALSE
  )
})

calculate_LCA <- function(contig_data){

  find_uniq <- function(x) {
    if (all(x=="")) return("")
    x_uniq <- unique(x[x!=""])
    result <- ifelse(length(x_uniq)==1,x_uniq,"-")
    return(result)
  }
  
  #Use find unique function to determine if certain taxon is "comon ancestor"
  contig_lca <- contig_data[,c("contig.id","file.id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom")][, lapply(.SD,find_uniq), .(contig.id, file.id)]
  
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
  
  #remove non LCA lineage
  contig_data[,c("tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"):=NULL]
  #Replace with LCA lineage
  contig_data <- merge.data.table(contig_data, contig_lca, by=c("contig.id","file.id"))
  
  return(contig_data)
}

annot.heatmap.data <- eventReactive({merged_annotationdata(); input$update_heatmap}, ignoreNULL = FALSE, {
  validate(
    need(!is.null(merged_annotationdata()),"Loading merged data for heatmap..."),
    need(!is.null(input$contig_table_rows_all)|length(input$contig_table_rows_all)>0,'Press the "Apply filters" button to load the heatmap')
  )
  
  merged_data <- merged_annotationdata()
  merged_data <- merged_data[input$contig_table_rows_all,]
  
  merged_data <- calculate_LCA(merged_data)
  
  #Remove columns not needed for plotting
  merged_data[,c("target","hit.start","hit.end","length.homology","identity","evalue(-log10)","bitscore","fraction.homology"):=NULL]
  #Unique rows to avoid counting contigs twice
  merged_data <- unique(merged_data)
  
  return(merged_data)
})

#Heatmap based on selected input
output$annot.heatmap <- renderRbokeh({
	
  req(annot.heatmap.data())
  
  x_axis <- req(input$metadata_stratify_selection)
  y_axis <- req(input$annotation_stratify_selection)
  
  fill <- req(input$annotation_count_selection)
  
  #Switch annotation based on slider
  switch (fill,
    countcontig = {
      #Sum contig by stratification
      plotdata <- annot.heatmap.data()[,.N, by=.(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    },
    abs_readcount = {
      #Sum reads per stratification
      plotdata <- annot.heatmap.data()[,sum(readcount), by=.(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    },
    scaled_readcount = {
      #Sum reads per stratification and scale by contig length
      plotdata <- annot.heatmap.data()[,sum(readcount/totalreadcount), by=.(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    }
  )
  
  if (input$logtransform) {
    #Log10 transform inplace
    plotdata[,(fill):=log10(.SD), .SDcols = fill]
  }
  
  x.axis <- annot.heatmap.data()[,unique(get(x_axis))]
	y.axis <- annot.heatmap.data()[,unique(get(y_axis))]
	
	colnames(plotdata)[colnames(plotdata)==fill] <- "Legend"
  
	plotdata <- plotdata[order(get(x_axis),get(y_axis))]
	
	figure(width = 600, height = 600,
	  xlab = x_axis, ylab = y_axis,
		tools = c("box_select", "save"),
		toolbar_location = "above") %>%
		ly_crect(data = plotdata, x=get(x_axis), y=get(y_axis), color=Legend, hover=list(x_axis=get(x_axis),y_axis=get(y_axis),fill=Legend), lname = "rects") %>%
		x_range(x.axis) %>%
		theme_axis("x", major_label_orientation = 45) %>%
		theme_grid("x", grid_line_alpha = 0) %>%
		y_range(y.axis) %>%
		theme_grid("y", grid_line_alpha = 0) %>%
		set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
		tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
})

merged_annotationdata <- reactive({
  validate(
    need(!is.null(annotation_data()),"Loading annotation data..."),
    need(!is.null(readcount_data()),"Loading readcount data..."),
    need(!is.null(metadata()|nrow(metadata()>0)),"Loading metadata...")
  )
  
  merged_data <- annotation_data()

  contig_length <- rbindlist(lapply(names(contig_data()), function(x){
    result <- seqinfo(contig_data()[[x]]$idx)
    result <- data.table("file.id"=x,"contig.id"=seqnames(result),"contig.length"=seqlengths(result))
  }))
  
  merged_data <- merge.data.table(merged_data,contig_length,by=c("file.id","contig.id"))
  
  merged_data[,fraction.homology:=length.homology/contig.length]
  
  if (isTruthy(readcount_data())) {
    merged_data <- merge.data.table(merged_data,readcount_data(),by=c("file.id","contig.id"))
  }
  
  meta_data <- metadata()
  
  #If selection is made based on metadata apply filter
  if (!is.null(input$metadata_filter_selection)) {
    meta_data <- metadata()[input$metadata_table_rows_all,]
  }
  
  merged_data <- merge.data.table(meta_data, merged_data, by='file.id')
  
  return(merged_data)
})

#Generate a table with contig info
output$contig_table <- renderDataTable({
  validate(
    need(!is.null(annotation_data()),"Loading annotation data..."),
    need(!is.null(merged_annotationdata()),"Loading merged data...")
  )
  
  contig_data <- merged_annotationdata()
  
  #Display only annotation quality filters
  filter <- colnames(contig_data)[colnames(contig_data)%in%c("readcount","contig.length","length.homology","fraction.homology","identity","evalue(-log10)","bitscore")]
  
  #Set the default quality threshold values
  searchCols <- rep(list(NULL), length(filter))
  names(searchCols) <- filter
  searchCols[names(default_values)] <- default_values
  names(searchCols) <- NULL
  searchCols <- c(list(NULL), searchCols)

  plotdata <- contig_data[,filter, with=F]
  datatable(plotdata,
            filter = 'top',
            options = list(
              searchCols = searchCols,
              dom = "t",
              pageLength = 0,
              # sScrollY = "0px",
              autoWidth = TRUE
            )
  )
})

#Generate a table with metadata info
output$metadata_table <- renderDataTable({
  validate(
    need(!is.null(metadata()),"Loading metadata...")
  )
  req(input$metadata_filter_selection)
  plotdata <- metadata()[,input$metadata_filter_selection, with=F]
  datatable(plotdata,
                filter = list(
                  position = 'top', 
                  clear = FALSE, 
                  plain = TRUE),
                options = list(
                  dom = "t",
                  pageLength = 0,
                  # sScrollY = "0px",
                  autoWidth = TRUE,
                  columnDefs = list(list(className = 'dt-left', targets = "_all"))
                  )
  )
})

#Selected fields from heatmap
heatmap.selected <- reactive({
  selection_data <- annot.heatmap.data()
  
  if (is.null(input$heatmap.select)) {
    selected_rows <- selection_data
    setkey(selected_rows, "contig.id")
  } else {
    x_axis <- req(input$metadata_stratify_selection)
    y_axis <- req(input$annotation_stratify_selection)
    
    selected <- selection_data[,.N, by=.(get(x_axis),get(y_axis))]
    names(selected) <- c(x_axis, y_axis, 'N')
    selected <- selected[order(get(x_axis),get(y_axis))]
    selected <- selected[input$heatmap.select+1]
    
    selected_rows <- selection_data[get(x_axis) %in% selected[[1]] & get(y_axis) %in% selected[[2]]]
    setkey(selected_rows, "contig.id")
  }
  return(selected_rows)
})