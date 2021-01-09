shinyServer(function(input, output, session) {
	
	importdata <- callModule(importFiles, "import")

	annotation_data <- importdata$annot
	default_values <- importdata$defaults
	contig_data <- importdata$contig
	readcount_data <- importdata$readcounts
  metadata <- importdata$metadata
  
	output$importdata <- renderMenu({
	  menuItem("File Import", tabName = "dataimport", selected=T)
	})
	
	source("heatmap_panel.R", local = TRUE)
	output$heatmap <- renderMenu({
		req(annotation_data(), contig_data(), readcount_data())
    
		menuItem("Interactive Data Browser", tabName = "heatmap")
	})
  
	source("seq_info.R", local = TRUE)
	output$seqinfo <- renderMenu({
		req(annotation_data(), contig_data())

		menuItem("Selected sequence annotations", tabName = "seqinfo")
	})
	
	isolate({updateTabItems(session, "sidebar", "dataimport")})
})
