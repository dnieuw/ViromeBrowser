#Import data tab#####
dataimport.tab <- tabItem(
  tabName = "dataimport",
  importFilesOutput("import")
)
#####

#Heatmap tab#####
heatmap.tab <- tabItem(
	tabName = "heatmap",
	fluidPage(
	  
	  fluidRow(
	    box(width=12, solidHeader = F, title = "Heatmap of lowest common ancestor annotations, based on selected filters", status = "primary",
	        fluidRow(
	          column(2,
	                 box(width=12,
	                     uiOutput("metadata_stratify_selection"),
	                     uiOutput("annotation_count_selection"),
	                     uiOutput("annotation_stratify_selection"),
	                     prettySwitch("logtransform","Log10 transform")
	                 )
	          ),
	          column(10, align="center",
	                 withSpinner(plotlyOutput("annot.heatmap", width = "50%", height = "800px")),
	                 actionButton("deselect_button","Deselect")
	          )
	        )
	    )
	  ),
	  fluidRow(
	    box(width = 12, solidHeader = T, title = "Metadata Filter Settings", 
	        status = "primary", collapsible = T, collapsed = F,
	        fluidPage(
	          DT::dataTableOutput("metadata_table")
	        )
	    )
	  ),
	  fluidRow(
	    box(width = 12, solidHeader = T, title = "Annotation Filter Settings",
	        collapsible = T, collapsed = F, status = "warning",
	          column(12, align="center",
               fluidPage(
                 DT::dataTableOutput("contig_table", width = "100%")
               )
	          )
	    )
	  )
	)
)
#####

#Seqinfo tab#####
seqinfo.tab <- tabItem(
	tabName = "seqinfo",
	fluidRow(
  	tabBox(width = 12,
  		tabPanel(title = "Annotation Table",
  		  fluidPage(
          column(2,
            box(width=NULL,
                h4("Applied filters:"),
                uiOutput("filterInfo"),
                br(),
                downloadButton("downloadContig", "Download Selected Contigs"),
                downloadButton("downloadAllContig", "Download All Filtered Contigs")
            )
          ),
          column(10,
            annotationTableOutput("annot.table")
          )
  		  )
  		),
  		###
  		tabPanel(title = "Contig Information",
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
  					    plotlyOutput("orfplot", width = "80%", height = "600px")
  					  )
  					)
  				)
  			)
  		),
  		###
  		tabPanel(title = "ORF Information",
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
  		tabPanel(title = "ORF Collection Table",
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
  a(href = "https://github.com/dnieuw/ViromeBrowser", target = "_blank", img(src="ViromeBrowser_logo.svg", align = "middle", width="100%")),
  sidebarMenu(id = "sidebar",
    menuItemOutput("importdata"),
    menuItemOutput("heatmap"),
    menuItemOutput("seqinfo")
  ),
  tags$footer(
    a(href = "https://www.compare-europe.eu/", target = "_blank", img(src="COMPARE_logo.png", width = "100%")),
    style=paste0("position:absolute; align: center; bottom:0px; width:100%; height:",
                 "110px", "; color: white; padding: 5px;")
  )
)



dashboardPage(skin = "green",
              dashboardHeader(title = "ViromeBrowser"),
              sidebar,
              body
)