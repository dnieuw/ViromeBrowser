
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
              "Absolute read count" = "abs_readcount",
              "Total scaled read count" = "scaled_readcount")

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
    x[x==""] <- "unassigned"
    if (all(x=="unassigned")) return("unassigned")
    x_uniq <- unique(x[x!="unassigned"])
    result <- ifelse(length(x_uniq)==1,x_uniq,"-")
    return(result)
  }
  
  #Use find unique function to determine if certain taxon is "comon ancestor"
  contig_lca <- contig_data[,c("contig.id","file.id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom")][, lapply(.SD,find_uniq), list(contig.id, file.id)]
  
  #Remove lineages with nonconforming lower ancestor
  contig_lca[superkingdom=='-',superkingdom:="unassigned"]
  contig_lca[superkingdom=='-',c("tax_name","species","genus","family","order","class","phylum","kingdom"):="unassigned"]
  contig_lca[kingdom=='-',c("tax_name","species","genus","family","order","class","phylum","kingdom"):="unassigned"]
  contig_lca[phylum=='-',c("tax_name","species","genus","family","order","class","phylum"):="unassigned"]
  contig_lca[class=='-',c("tax_name","species","genus","family","order","class"):="unassigned"]
  contig_lca[order=='-',c("tax_name","species","genus","family","order"):="unassigned"]
  contig_lca[family=='-',c("tax_name","species","genus","family"):="unassigned"]
  contig_lca[genus=='-',c("tax_name","species","genus"):="unassigned"]
  contig_lca[species=='-',c("tax_name","species"):="unassigned"]
  contig_lca[tax_name=='-',c("tax_name"):="unassigned"]
  
  #remove non LCA lineage
  contig_data[,c("tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"):=NULL]
  #Replace with LCA lineage
  contig_data <- merge.data.table(contig_data, contig_lca, by=c("contig.id","file.id"))
  
  return(contig_data)
}

annot.heatmap.data <- reactive({
  validate(
    need(!is.null(selected_data()),"Loading selected data..")
  )

  merged_data <- calculate_LCA(selected_data())
  
  #Remove columns not needed for plotting
  merged_data[,c("target","hit.start","hit.end","length.homology","identity","evalue(-log10)","bitscore","fraction.homology"):=NULL]
  #Unique rows to avoid counting contigs twice
  merged_data <- unique(merged_data)
  
  if(nrow(merged_data)==0) return(NULL)
  
  return(merged_data)
})

annot.heatmap.data.debounced <- annot.heatmap.data %>% debounce(1000)

#Heatmap based on selected input
output$annot.heatmap <- renderRbokeh({

  validate(
    need(!is.null(annot.heatmap.data.debounced()),"No data for currently selected filters..")
  )
  
  x_axis <- req(input$metadata_stratify_selection)
  y_axis <- req(input$annotation_stratify_selection)
  fill <- req(input$annotation_count_selection)

  #Switch annotation based on slider
  switch (fill,
    countcontig = {
      #Sum contig by stratification
      plotdata <- annot.heatmap.data.debounced()[,.N, by=list(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    },
    abs_readcount = {
      #Sum reads per stratification
      plotdata <- annot.heatmap.data.debounced()[,sum(readcount), by=list(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    },
    scaled_readcount = {
      #Sum reads per stratification and scale by contig length
      plotdata <- annot.heatmap.data.debounced()[,sum(readcount/totalreadcount), by=list(get(x_axis),get(y_axis))]
      names(plotdata) <- c(x_axis, y_axis, fill)
    }
  )
  
  if (input$logtransform) {
    #Log10 transform inplace
    plotdata[,(fill):=log10(.SD), .SDcols = fill]
  }
  
  x.axis <- annot.heatmap.data.debounced()[,unique(get(x_axis))]
	y.axis <- annot.heatmap.data.debounced()[,unique(get(y_axis))]
	
	plotdata <- plotdata[order(get(x_axis),get(y_axis))]
	
	switch (fill,
	        countcontig = {
	          hover <- list()
	          colnames(plotdata)[colnames(plotdata)==fill] <- "Contig count"
	          figure(width = 600, height = 600,
	                 xlab = x_axis, ylab = y_axis,
	                 tools = c("box_select", "save"),
	                 toolbar_location = "above") %>%
	            ly_crect(data = plotdata, x=get(x_axis), y=get(y_axis), color=`Contig count`, hover=list(x_axis=get(x_axis),y_axis=get(y_axis),'Contig count'=`Contig count`), lname = "rects") %>%
	            x_range(x.axis) %>%
	            theme_axis("x", major_label_orientation = 45) %>%
	            theme_grid("x", grid_line_alpha = 0) %>%
	            y_range(y.axis) %>%
	            theme_grid("y", grid_line_alpha = 0) %>%
	            set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
	            tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
	         },
	        abs_readcount = {
	          colnames(plotdata)[colnames(plotdata)==fill] <- "Absolute read count"
	          figure(width = 600, height = 600,
	                 xlab = x_axis, ylab = y_axis,
	                 tools = c("box_select", "save"),
	                 toolbar_location = "above") %>%
	            ly_crect(data = plotdata, x=get(x_axis), y=get(y_axis), color=`Absolute read count`, hover=list(x_axis=get(x_axis),y_axis=get(y_axis),'Absolute read count'=`Absolute read count`), lname = "rects") %>%
	            x_range(x.axis) %>%
	            theme_axis("x", major_label_orientation = 45) %>%
	            theme_grid("x", grid_line_alpha = 0) %>%
	            y_range(y.axis) %>%
	            theme_grid("y", grid_line_alpha = 0) %>%
	            set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
	            tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
	          },
	        scaled_readcount = {
	          colnames(plotdata)[colnames(plotdata)==fill] <- "Total scaled read count"
	          figure(width = 600, height = 600,
	                 xlab = x_axis, ylab = y_axis,
	                 tools = c("box_select", "save"),
	                 toolbar_location = "above") %>%
	            ly_crect(data = plotdata, x=get(x_axis), y=get(y_axis), color=`Total scaled read count`, hover=list(x_axis=get(x_axis),y_axis=get(y_axis),'Total scaled read count'=`Total scaled read count`), lname = "rects") %>%
	            x_range(x.axis) %>%
	            theme_axis("x", major_label_orientation = 45) %>%
	            theme_grid("x", grid_line_alpha = 0) %>%
	            y_range(y.axis) %>%
	            theme_grid("y", grid_line_alpha = 0) %>%
	            set_palette(continuous_color = pal_gradient(cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))), pal_size(min=0,max=2500)) %>%
	            tool_box_select(callback = shiny_callback("heatmap.select"),ref_layer = "rects")
	          }
	)
	
	
})

merged_annotationdata <- reactive({
  validate(
    need(!is.null(annotation_data()),"Loading annotation data..."),
    need(!is.null(readcount_data()),"Loading readcount data...")
  )
  
  merged_data <- annotation_data()

  contig_length <- rbindlist(lapply(names(contig_data()), function(x){
    result <- seqinfo(contig_data()[[x]]$idx)
    result <- data.table("file.id"=x,"contig.id"=seqnames(result),"contig.length"=seqlengths(result))
  }))
  
  merged_data <- merge.data.table(merged_data,contig_length,by=c("file.id","contig.id"))
  merged_data[,fraction.homology:=length.homology/contig.length]
  
  merged_data <- merge.data.table(merged_data,readcount_data(),by=c("file.id","contig.id"))
  
  if(nrow(merged_data)==0){
    return(NULL)
  }
  
  return(merged_data)
})

selected_annotationdata <- reactive({
  validate(
    need(!is.null(merged_annotationdata()),"Loading annotation data...")
  )
  annotationdata <- merged_annotationdata()[input$contig_table_rows_all]
  
  if (nrow(annotationdata)==0) {
    return(NULL)
  }
  
  return(annotationdata)
})

selected_metadata <- reactive({
  validate(
    need(!is.null(metadata()),"Loading metadata...")
  )
  metadata <- metadata()[input$metadata_table_rows_all]
  
  if (nrow(metadata)==0) {
    return(NULL)
  }
  
  return(metadata)
})

selected_data <- reactive({
  validate(
    need(!is.null(selected_metadata()),"No contigs for selected metadata..."),
    need(!is.null(selected_annotationdata()),"No contigs for selected contig filters...")
  )
  return(merge.data.table(selected_metadata(),selected_annotationdata(), by = "file.id"))
})

#Generate a table with contig info
output$contig_table <- renderDataTable({
  validate(
    need(!is.null(merged_annotationdata()),"Loading annotation data...")
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
  activate <- input$annotation_stratify_selection#No idea why, but this is needed to have the heatmap show up right away
  plotdata <- contig_data[,filter, with=F]
  rownames(plotdata) <- contig_data[,list(unique_names=paste(contig.id,seq(.N),sep="#")),contig.id][,unique_names]
  datatable(plotdata,
            rownames = T,
            filter = 'top',
            options = list(
              searchCols = searchCols,
              ordering=T,
              dom = "tip",
              pageLength = 10,
              scrollX = TRUE
            )
  )
})

#Generate a table with metadata info
output$metadata_table <- renderDataTable({
  validate(
    need(!is.null(metadata()),"Loading metadata...")
  )
  plotdata <- metadata()
  datatable(plotdata,
                filter = list(
                  position = 'top', 
                  clear = FALSE, 
                  plain = TRUE),
                  options = list(
                    ordering=T,
                    dom = "tip",
                    pageLength = 10,
                    scrollX = TRUE,
                    scrollY = TRUE,
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
    
    selected <- selection_data[,.N, by=list(get(x_axis),get(y_axis))]
    names(selected) <- c(x_axis, y_axis, 'N')
    selected <- selected[order(get(x_axis),get(y_axis))]
    selected <- selected[input$heatmap.select+1]
    
    selected_rows <- selection_data[get(x_axis) %in% selected[[1]] & get(y_axis) %in% selected[[2]]]
    setkey(selected_rows, "contig.id")
  }
  return(selected_rows)
})