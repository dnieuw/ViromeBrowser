shinyServer(function(input, output, session) {
	
	importdata <- callModule(importFiles, "import")

	annotation_data <- importdata$annot
	default_values <- importdata$defaults
	contig_data <- importdata$contig
  metadata <- importdata$metadata
  
	output$importdata <- renderMenu({
	  menuItem("File Import", tabName = "dataimport", selected=T)
	})
	
	source("heatmap_panel.R", local = TRUE)
	output$heatmap <- renderMenu({
		req(annotation_data())

		menuItem("Interactive Data Browser", tabName = "heatmap")
	})
  
	source("seq_info.R", local = TRUE)
	output$seqinfo <- renderMenu({
		req(annotation_data(), contig_data())

		menuItem("Sequence Information", tabName = "seqinfo")
	})
	
	isolate({updateTabItems(session, "sidebar", "dataimport")})
})
