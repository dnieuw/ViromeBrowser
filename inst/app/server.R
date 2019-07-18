shinyServer(function(input, output, session) {
	
	importdata <- callModule(importFiles, "EBI")

	annotation_data <- importdata$annot
	contig_data <- importdata$contig

	output$importdata <- renderMenu({
	  menuItem("Import data", tabName = "dataimport", selected=T)
	})
	
	source("heatmap_panel.R", local = TRUE)
	output$heatmap <- renderMenu({
		req(annotation_data())

		menuItem("Sample heatmap", tabName = "heatmap")
	})
  
	source("seq_info.R", local = TRUE)
	output$seqinfo <- renderMenu({
		req(annotation_data(), contig_data())

		menuItem("Sequence information", tabName = "seqinfo")
	})
	
	isolate({updateTabItems(session, "sidebar", "dataimport")})
})
