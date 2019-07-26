#' @title Launch the Virome Browser
#'
#' @description
#' An application for viewing virome annotation files
#'
#' @import shiny
#' @import ggplot2
#' @import shinydashboard
#' @import DT
#' @import rbokeh
#' @import Biostrings
#' @import Rsamtools
#' @import data.table
#' @import markdown
#' @import shinyWidgets
#' @import stringr
#' @export
#' @param host The host ip address, 127.0.0.1 or "localhost" by default.
#' @param port The host port, 3838 by default.
#' @examples 
#' \donttest{
#' viromeBrowser()
#' }
viromeBrowser <- function(host = getOption("shiny.host", "127.0.0.1"), port = getOption("shiny.port",3838)) {
	appDir = system.file('app', package='viromeBrowser', mustWork=TRUE)
	runApp(appDir = appDir,
		host = host,
		port = port,
		workerId = "",
		quiet = FALSE,
		display.mode = "normal",
		test.mode = FALSE)
}
