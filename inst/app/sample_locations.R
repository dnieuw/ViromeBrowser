#Rscript to generate a leaflet widget of sample locations

output$graph.locations <- renderLeaflet({
	metadata <- metadata$imported()
	leaflet(metadata) %>%
	addTiles() %>%
	addCircles(lng= ~ longitude, 
		lat= ~ latitude, 
		radius = ~demography/10, 
		layerId = ~ city)
})


# Show a popup at the given location
showPopup <- function(id, lat, lng) {
	metadata <- metadata$imported()
  content <- as.character(tagList(
    tags$h4(metadata[city==id,city]),
    tags$strong(metadata[city==id,site]),
  	tags$br(),
  	tags$strong(metadata[city==id,site.characteristic]),
  	tags$br(),
    sprintf("Estimated coverage: %s km^2", metadata[city==id,site.size]),
  	tags$br(),
    sprintf("Estimated number of people served: %s million", metadata[city==id,demography]/1000000)
  ))
  leafletProxy("graph.locations") %>% addPopups(lng, lat, content, layerId = id)
}

# Register marker click event and show popup accordingly
observe({
  leafletProxy("graph.locations") %>% clearPopups()
  event <- input$graph.locations_shape_click
  if (is.null(event))
    return()

  isolate({
    showPopup(event$id, event$lat, event$lng)
  })
})

# Remove popups when clicked on map
observe({
  leafletProxy("graph.locations") %>% clearPopups()
  event <- input$graph.locations_click
})


output$selected.city <- renderPrint({
  event <- input$graph.locations_shape_click
  if (is.null(event))
    return()
  return(event)
})