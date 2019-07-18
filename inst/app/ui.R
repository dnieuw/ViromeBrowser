#Import data tab#####
dataimport.tab <- tabItem(
  tabName = "dataimport",
  importFilesOutput("EBI")
)
#####

#Heatmap tab#####
heatmap.tab <- tabItem(
	tabName = "heatmap",
	fluidPage(
		fluidRow(
			box(width = 12,
				uiOutput("select.tax")
			)
		),
		fluidRow(
		  box(width = 12,
		    DT::dataTableOutput("table")
		  )
		),
		fluidRow(
		  box(width=12,
  		  column(12, align="center",
          rbokehOutput("annot.heatmap", width = "50%", height = "600px")
  		  )
  		)
		)
	)
)
#####

#Annotation table tab#####
annotationTableOutput <- function(id) {
	ns <- NS(id)
	fluidPage(
		box(width = 12,
			DT::dataTableOutput(ns("table")),
			actionButton(ns("clearbutton"), "Clear selection")
		)
	)
}
#####

#Seqinfo tab#####
seqinfo.tab <- tabItem(
	tabName = "seqinfo",
	fluidRow(
  	tabBox(width = 12,
  		tabPanel(title = "Annotation table",
  			annotationTableOutput("annot.table"),
  			downloadButton("downloadContig", "Download Selected")
  		),
  		###
  		tabPanel(title = "Contig information",
  			fluidPage(
  				fluidRow(
  					box(status = "primary",
  						width = 12,
  						solidHeader = T,
  						fluidRow(
  							column(6,
  								uiOutput("select.contig.current"),
  								uiOutput("contig.summary")
  							),
  							column(6,
  								uiOutput("orf.size.cutoff")
  							)
  						)
  					)
  				),
  				fluidRow(
  					box(width = 12,
  					  column(12, align="center",
  					    rbokehOutput("orfplot", width = "80%", height = "600px")
  					  )
  					)
  				)
  			)
  		),
  		###
  		tabPanel(title = "ORF information",
  			fluidPage(
  				fluidRow(
  					box(status = "primary",
  						width = 12,
  						solidHeader = T,
  						fluidRow(
  							column(6,
  								uiOutput("select.orf.current"),
  								uiOutput("orf.sequence"),
  								uiOutput("orf.sequence.type"),
  								actionButton("collect.orf","Collect ORFs")
  							),
  							column(6,
  								uiOutput("orf.slidingwindow.winsize"),
  								uiOutput("orf.slidingwindow.type"),
  								tags$img(src="http://www.jalview.org/help/html/misc/properties.gif", height = "250px")
  							)
  						)
  					)
  				),
  				fluidRow(
  					box(width = 12,
  						plotOutput("orf.slidingwindow")
  					)
  				)
  			)
  		),
  		###
  		tabPanel(title = "Sequence collection table",
  			fluidPage(
  				fluidRow(
  					box(status = "primary",
  						width = 12,
  						solidHeader = T,
  						DT::dataTableOutput("orf.collection.table"),
  						uiOutput("download.seqtype"),
  						downloadButton('downloadFasta', 'Download Fasta'),
  						actionButton("clear.orf","Clear Collection")
  					)
  				)
  			)
  		)
  	)
	)
)
#####

body <- dashboardBody(
  tags$script(
    "Shiny.addCustomMessageHandler(
			'resetValue', function(variableName) {
				Shiny.onInputChange(variableName, null);
			}
		);"
  ),
  tabItems(
    dataimport.tab,
    heatmap.tab,
    seqinfo.tab
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(id = "sidebar",
    menuItemOutput("importdata"),
    menuItemOutput("heatmap"),
    menuItemOutput("seqinfo")
  )
)

dashboardPage(skin = "green",
              dashboardHeader(title = "ViromeBrowser"),
              sidebar,
              body
)