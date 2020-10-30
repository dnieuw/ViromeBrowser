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

output$annotation_stratify_selection <- renderUI({
  
  annotation_levels <- c("superkingdom",
                         "kingdom",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "species",
                         "taxonomic_name")
  
  available_annotations <- annotation_levels[annotation_levels %in% colnames(annotation_data())]
  
  pickerInput(
    inputId = "annotation_stratify_selection",
    label = "Annotate with:", 
    choices = available_annotations,
    options =  list(
      size = 5
    ), 
    multiple = FALSE
  )
})

annot.heatmap.data <- eventReactive({annotation_data(); input$update_heatmap}, ignoreNULL = FALSE, {
  validate(
    need(!is.null(annotation_data()),"Loading..."),
    need(!is.null(metadata()),"Loading...")
  )
  
  contig_data <- annotation_data()
  
  #If no selection is made based on contig data set the default thresholds
  if (input$contig_table_set_defaults) {
    
    #Initialize default values
    sapply(names(default_values), function(name){
      default <- default_values[[name]]

      if (is.na(default[1])){
        #Skip non-set defaults
        return(NULL)
      }
      if(is.numeric(default)){
        #It is a slider range
        contig_data <<- contig_data[get(name)>default[1]&get(name)<default[2]]
        return(NULL)
      }
      #It is a selection of factors/characters
      contig_data <<- contig_data[get(name)%in%default]
      return(NULL)
    })
  } else {
    contig_data <- contig_data[input$contig_table_rows_all,]
  }
  
  meta_data <- metadata()
  
  #If selection is made based on metadata apply filter
  if (!is.null(input$metadata_filter_selection)) {
    meta_data <- metadata()[input$metadata_table_rows_all,]
  }
  
  merged_data <- merge(contig_data,meta_data, by='file.id')

  return(merged_data)
})

#Heatmap based on selected input
output$annot.heatmap <- renderRbokeh({
	
  x_axis <- req(input$metadata_stratify_selection)
  y_axis <- req(input$annotation_stratify_selection)
  
  #Sum contig by stratification
	plotdata <- annot.heatmap.data()[,.N, by=.(get(x_axis),get(y_axis))]
	names(plotdata) <- c(x_axis, y_axis, 'N')
	
	x.axis <- annot.heatmap.data()[,unique(get(x_axis))]
	y.axis <- annot.heatmap.data()[,unique(get(y_axis))]
	
	Contigs <- plotdata[,N]
	
	plotdata <- plotdata[order(get(x_axis),get(y_axis))]
	
	figure(width = 600, height = 600,
	  xlab = x_axis, ylab = y_axis,
		legend_location = NULL,
		tools = c("box_select", "save"),
		toolbar_location = "above") %>%
		theme_legend(background_fill_alpha = 0, border_line_alpha = 0) %>%
		ly_crect(data = plotdata, x=get(x_axis), y=get(y_axis), color=N, hover=list(x_axis=get(x_axis),y_axis=get(y_axis),"#Contigs"=N), lname = "rects") %>%
		x_range(x.axis) %>%
		theme_axis("x", major_label_orientation = 45) %>%
		theme_grid("x", grid_line_alpha = 0) %>%
		y_range(y.axis) %>%
		theme_grid("y", grid_line_alpha = 0) %>%
		set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
		tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
})

#Generate a table with contig info
output$contig_table <- renderDataTable({
  validate(
    need(!is.null(annotation_data()),"Loading...")
  )
  plotdata <- annotation_data()[,c(-1,-2), with=FALSE]
  datatable(plotdata,
            filter = 'top',
            options = list(
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
    need(!is.null(metadata()),"Loading...")
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