output$select.tax <- renderUI({
		opts <- req(annotation_data())[,unique(annotation)]
		selectizeInput("select.tax", "Select taxonomy to display:", opts,
		selected = opts[1:10],
		width= "100%",
		multiple=TRUE)
})

annot.tax <- reactive({
	req(input$select.tax)
	annot <- annotation_data()[annotation%in%input$select.tax]
	return(annot)
})

annot.heatmap.data <- eventReactive(input$table_rows_all,{
  validate(
    need(!is.null(annot.tax()),"Loading..."),
    need(!is.null(input$table_rows_all), "Loading...")
  )
  selected_rows <- input$table_rows_all
  
  annot_sub <- annot.tax()[selected_rows,]
  
  return(annot_sub)
})

#Heatmap based on selected input
output$annot.heatmap <- renderRbokeh({
	#Sum contig by file.id and annotation
	plotdata <- req(annot.heatmap.data())[,.(N=length(unique(contig.id))),by=.(file.id,annotation)]
	
	x.axis <- annotation_data()[,unique(file.id)]
	y.axis <- input$select.tax
	Contigs <- plotdata[,N]
	figure(width = 600, height = 600,
	  xlab = "Filenames", ylab = "Annotation",
		legend_location = NULL,
		tools = c("box_select"),
		toolbar_location = "above") %>%
		theme_legend(background_fill_alpha = 0, border_line_alpha = 0) %>%
		ly_crect(data = plotdata, x=file.id, y=annotation, color=Contigs, hover=list("File"=file.id,"Annotation"=annotation,"#Contigs"=N), lname = "rects") %>%
		x_range(x.axis) %>%
		theme_axis("x", major_label_orientation = 45) %>%
		theme_grid("x", grid_line_alpha = 0) %>%
		y_range(y.axis) %>%
		theme_grid("y", grid_line_alpha = 0) %>%
		set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
		tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
})

heatmap.selected <- reactive({
  selected <- isolate(annot.heatmap.data())
  selected <- selected[,.(N=length(unique(contig.id))),by=.(file.id,annotation)]
	if (is.null(input$heatmap.select))
		selected <- selected
	else
		selected <- selected[input$heatmap.select+1]
	return(selected)
})

#Generate a table with contig info
output$table <- DT::renderDataTable({
  plotdata <- req(annot.tax()[,c("contig.length","length.homology","frac.homology","identity","evalue")])
  DT::datatable(plotdata,
            filter = 'top',
            selection = list(
              mode = 'multiple',
              selected = 1,
              target = 'row'
            ),
            options = list(
              ordering=F,
              dom = "ti",
              pageLength = 0,
              sScrollY = "0px",
              autoWidth = TRUE
            )
  )
})
