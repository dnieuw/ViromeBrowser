find.ORFs <- function(contig, stopcodon = c("TAA","TGA","TAG")){
	stopcodon <- DNAStringSet(stopcodon)
	#Match stopcodon to sequence in forward direction
	matches.fwd <- do.call(c,lapply(as.character(stopcodon),matchPattern,contig))
	mcols(matches.fwd) <- data.frame("frame"=paste0("+",end(matches.fwd)%%3+1))

	#Match stopcodon to sequence in reverse direction
	matches.rev <- do.call(c,lapply(as.character(reverseComplement(stopcodon)),matchPattern,contig))
	mcols(matches.rev) <- data.frame("frame"=paste0("-",end(matches.rev)%%3+1))

	matches <- c(matches.fwd,matches.rev)

		allorfs <- do.call(c,lapply(paste0(c("+","-"),c(1,1,2,2,3,3)), function(x){
		#Find gaps between fwd stopcodons in each frame
		orfs <- gaps(matches[mcols(matches)$frame==x])
		mcols(orfs)$frame <- as.character(x)
		#Trim 1 and 2 nt ORFS
		orfs <- orfs[width(orfs)>2]
		#Continue if no ORFs are found
		if(length(orfs)==0) {
		  return(NULL)
		}
		#Trim start and end for residual nt
		start(orfs[1]) <- start(orfs[1])+width(orfs[1])%%3
		end(orfs[length(orfs)]) <- end(orfs[length(orfs)])-width(orfs[length(orfs)])%%3
		return(orfs)
	}))
	mcols(allorfs)$orfid <- paste("ORF",seq(length(allorfs)),sep="")
	return(allorfs)
}

plot.orfs <- function(orfs,lim){
	plotdata <- as.data.table(ranges(orfs))
	limits.x <- c(-5000,max(end(orfs))+5000)
	limits.y <- c(-2.5,2.5)
	plotdata <- orfs[width(orfs)>=lim]
	Frame <- mcols(plotdata)$frame
	p <- figure(width = 800, height = 400, tools = c("pan","wheel_zoom","reset"), toolbar_location = "above", lod_threshold = 1000, xlim=limits.x, ylim=limits.y) %>%
	ly_annular_wedge(data = as.data.table(ranges(plotdata)),
		x = start+width/2, y = 0, inner_radius = width/2, outer_radius = width/1.95,
		start_angle = sapply(as.character(mcols(plotdata)$frame), function(x)switch(x,"-1"=pi,"-2"=pi,"-3"=pi,"+1"=0,"+2"=0,"+3"=0)),
		end_angle = sapply(as.character(mcols(plotdata)$frame), function(x)switch(x,"-1"=0,"-2"=0,"-3"=0,"+1"=pi,"+2"=pi,"+3"=pi)),
		direction = "anticlock", color = Frame, alpha = 1, lname = "wedges") %>%
		set_palette(discrete_color = pal_color(c(	"#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))) %>%
		y_axis(visible = F) %>% x_axis(label = "Stopcodon Position") %>%
		tool_wheel_zoom("width") %>%
		tool_pan("width") %>%
		tool_box_zoom()
	return(p)
}

plot.AAfreq <- function(sequence) {
	freq <- letterFrequency(sequence,AA_STANDARD)
	names(freq) <- AA_DATA[names(freq),triple]
	plotdata <- data.frame("amino"=factor(names(freq),names(freq),ordered = T) , "freq"=freq)

	ggplot(plotdata) +
		geom_polygon(aes(x=amino, y=freq,group=1), color = "purple", fill=NA) + coord_polar()
}

plot.slidingwindow <- function(sequence,winsize,type){
	alph <- switch(class(sequence),"DNAString"=DNA_BASES,"AAString"=AA_STANDARD)
	slide <- as.data.table(letterFrequencyInSlidingView(sequence,winsize,alph,as.prob = T))
	switch(type,
		"seq"={
			if(all(alph==AA_STANDARD))
				names(slide) <- AA_DATA[names(slide),full]
			slide[,"location":=1:nrow(slide)]
			slide.melt <- melt(slide, id.vars = "location", measure.vars = seq(length(alph)))
			ggplot(slide.melt) +
				geom_tile(aes(x=location,y=variable,fill=value)) +
				scale_fill_gradientn(colors=c("gray",rev(RColorBrewer::brewer.pal(11, "Spectral")))) +
				labs(x="Sliding window location", y="", fill="Relative\nabundance") +
				theme(panel.background = element_blank(),
					axis.text.y=element_text(size=16),
					axis.ticks.y=element_blank())
		},
		"prop"={
			plotdata <- data.table("location"=1:nrow(slide))
			plotdata[,names(AA_DATA)[4:18] := lapply(names(AA_DATA)[4:18], function(prop) rowSums(slide[,lapply(names(slide),function(x) get(x)*AA_DATA[x, get(prop)])]))]
			plotdata.melt <- melt(plotdata, id.vars = "location", measure.vars = 2:15)
			ggplot(plotdata.melt) +
				geom_tile(aes(x=location,y=variable,fill=value)) +
				scale_fill_gradientn(colors=c("gray",rev(RColorBrewer::brewer.pal(11, "Spectral")))) +
				labs(x="Sliding window location", y="", fill="Value") +
				theme(panel.background = element_blank(),
					axis.text.y=element_text(size=16),
					axis.ticks.y=element_blank())
		}
	)
}

annotationTable <- function(input, output, session, datasubset, alldata) {
  ns <- session$ns
  
	table.data <- reactive({

		selected <- datasubset()
		validate(need(!is.null(alldata()),message = "Loading..."))
		annot <- alldata()[annotation%in%selected[,annotation] & file.id%in%selected[,file.id]]
		setkey(annot,contig.id)
		return(annot)
	})
	
	#Generate a table with contig info
	output$table <- DT::renderDataTable({
		plotdata <- table.data()
		DT::datatable(plotdata,
			selection = list(
	  		mode = 'multiple',
	  		selected = 1,
	  		target = 'row'
	  	),
	  	options = list(
	  		ordering=F,
	  		dom = "tip",
	  		pageLength = 100,
	  		scrollX = TRUE,
	  		scrollY = 100,
	  		sScrollY = "500px"
	    )
	  )
	})
  
	proxy = DT::dataTableProxy('table')
	
	observeEvent(input$clearbutton, {
	  proxy %>% DT::selectRows(NULL)
	})
	
	get.selected <- reactive({
		return(input$table_rows_selected)
	})

	return(list("table" = table.data, "selected" = get.selected))
}

annot <- callModule(annotationTable, "annot.table", heatmap.selected, annot.heatmap.data)

get.contig.selected <- reactive({
	selected <- annot$table()[annot$selected()]
	validate(need(!nrow(selected)==0,message="Please select a contig to visualize"))
	return(selected)
})

output$select.contig.current <- renderUI({
	selected <- get.contig.selected()
	choices <- seq(selected[,contig.id])
	names(choices) <- selected[,contig.id]

	ui <- pickerInput(inputId = "current_contig",
  	label = "Select Contig",
		choices = choices,
		choicesOpt = list(subtext = paste("Annotation", selected[,annotation], sep = ": "))
	)
	return(ui)
})

get.contig.current <- reactive({
	selected <- get.contig.selected()
	current <- selected[as.numeric(req(input$current_contig))]
	return(current)
})

get.contig.seq <- reactive({
	current <- get.contig.current()
	contigid <- current$contig.id
	file <- current$file.id
  
	if (any(file == names(contig_data()))==FALSE) {
	  showModal(modalDialog(title = "ERROR!",
	              paste0("Could not find a fasta file for: ",file),
	              easyClose = T,
	              footer = NULL))
	  return(NULL)
	}

	fasta <- contig_data()[[file]]$fasta
	idx <- contig_data()[[file]]$idx
	

	if (any(seqlevels(idx)==contigid)==FALSE) {
	  showModal(modalDialog(title = "ERROR!",
	              paste0("Could not find contig: ",contigid," in your fasta file"),
	              easyClose = T,
	              footer = NULL))
	  return(NULL)
	}
	
	loc <- idx[seqlevels(idx)==contigid]
	
	contig <- scanFa(fasta, param=loc)[[1]]

	mcols(contig)$file <- file
	mcols(contig)$contigid <- contigid

	return(contig)
})

get.contig.orfs <- reactive({
  req(get.contig.seq())
	contig.seq <- get.contig.seq()
	orfs <- find.ORFs(contig.seq)
	mcols(orfs)$contigid <- mcols(contig.seq)$contigid[1]
	mcols(orfs)$file <- mcols(contig.seq)$file[1]
	return(orfs)
})

get.orf.selected <- reactive({
	orfs <- get.contig.orfs()
	orfid <- which(width(orfs)>=input$orf.size.cutoff)
	return(orfs[orfid])
})

get.current.orf <- reactive({
	req(input$orfs)
	orf <- get.orf.selected()
	orf <- orf[which(mcols(orf)$orfid==input$orfs)]
	return(orf)
})

get.orf.seq <- reactive({
	orf <- get.current.orf()

	orf.seq <- as(orf[[1]], "DNAString")

	if(mcols(orf)$frame%in%c("-1","-2","-3")){
		orf.seq <- reverseComplement(orf.seq)
	}

	if(input$orf.sequence.type){
		orf.seq <- translate(orf.seq)
	}

	return(orf.seq)
})

output$contig.summary <- renderUI({
	selected <- get.contig.current()
	fields <- tags$div(lapply(names(selected), function(field){
		list(tags$strong(paste(field,": ")), selected[,get(field)], br())
	}))
	return(fields)
})

output$orfplot <- renderRbokeh({
	selected <- get.contig.current()
	hit.start <- selected$hit.start
	hit.end <- selected$hit.end

	orfs <- get.contig.orfs()
	if (any(is.na(c(hit.start,hit.end))))
		plot.orfs(orfs,input$orf.size.cutoff)
	else
		plot.orfs(orfs,input$orf.size.cutoff) %>% ly_abline(v=hit.start) %>% ly_abline(v=hit.end)
})

output$orf.size.cutoff <- renderUI({
	orfs <- get.contig.orfs()
	min.orf.size <- min(width(orfs))
	max.orf.size <- max(width(orfs))
	sliderInput("orf.size.cutoff","Minimal Open Reading Frame size:", min = min.orf.size, max = max.orf.size, value = 250)
})

output$select.orf.current <- renderUI({
	orfs <- get.orf.selected()

	ranges <- as.data.frame(ranges(orfs))
	meta <- as.data.frame(mcols(orfs))
	data <- cbind(meta, ranges)

	x <- pickerInput(inputId = "orfs",
  	label = "Select ORF",
		choices = data$orfid,
  	choicesOpt = list(subtext = paste(
  		paste("frame", data$frame, sep = ": "),
  		paste("start", data$start, sep = ": "),
  		paste("end", data$end, sep = ": "),
  		paste("width", data$width, sep = ": ")
  	))
	)
	return(x)
})

output$orf.sequence.type <- renderUI({
	checkboxInput("orf.sequence.type", "Translate", value=F)
})

output$orf.sequence <- renderUI({
	return(tags$textarea(get.orf.seq(), rows="10", style="width:100%"))
})

output$orf.slidingwindow.winsize <- renderUI({
	max <- length(get.orf.seq())

	sliderInput("orf.slidingwindow.winsize", "Sliding window size:", min = 1, max=max, value = 50)
})

output$orf.slidingwindow.type <- renderUI({
	checkboxInput("orf.slidingwindow.type", "Properties", value=F)
})

output$orf.slidingwindow <- renderPlot({
	size <- input$orf.slidingwindow.winsize
	type <- ifelse(input$orf.slidingwindow.type,"prop","seq")

	validate(need(!is.null(size),message="Loading..."))

	if (input$orf.sequence.type){
		return(plot.slidingwindow(get.orf.seq(),size,type))
	} else {
		return(plot.slidingwindow(get.orf.seq(),size,"seq"))
	}
})

orf.collection <- reactiveValues(orfs = NULL)

observeEvent(input$collect.orf,{
	#Get all selected ORF from the plot
	selected.orf <- isolate(get.orf.selected())
	#Iteratively add them to the ORF collection
	lapply(seq_along(selected.orf), function(i){
		orf.collection$orfs <- c(orf.collection$orfs,selected.orf[i])
	})
})

observeEvent(input$clear.orf,{
		orf.collection$orfs <- NULL
})

output$orf.collection.table <-  DT::renderDataTable({
	ranges <- do.call(c,lapply(orf.collection$orfs,ranges))
	meta <- do.call(rbind,lapply(orf.collection$orfs,mcols))
	data <- cbind(as.data.frame(meta),as.data.frame(ranges))
	outputtable <- DT::datatable(data,
		selection = list(
  		mode = 'multiple',
  		selected = NULL,
  		target = 'row'
  	),
  	options = list(
  		ordering=F,
  		dom = "t",
  		pageLength = 100,
  		lengthMenu = 25,
  		scrollX = TRUE,
  		scrollY = 100,
  		sScrollY = "400px"
    )
  )
	return(outputtable)
})

output$download.seqtype <- renderUI({
	radioGroupButtons(inputId = "download.seqtype", label = "Sequence type", choices = c("Nucleotide", "Aminoacid"),
		checkIcon = list(
			yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
    	no = tags$i(class = "fa fa-square-o", style = "color: steelblue")
		)
	)
	# radioGroupButtons(inputId = "download.seqtype", label = "Sequence type", choices = c("Nucleotide", "Aminoacid"), justified = TRUE)
})

output$downloadFasta <- downloadHandler(
	filename = function() {
		filename = paste("fasta", ".fa", sep='')
	},
	content = function(file) {
		lapply(orf.collection$orfs, function(orf) {
			sequence <- DNAStringSet(orf)
			if(mcols(orf)$frame%in%c("-1","-2","-3")){
				sequence <- reverseComplement(sequence)
			}
			if(input$download.seqtype=="Aminoacid"){
				sequence <- translate(sequence)
			}
			names(sequence) <- paste(unlist(mcols(orf)), collapse = "_")
			writeXStringSet(sequence,file,append=T)
		})
	}
)

output$downloadContig <- downloadHandler(
	filename = function() {
		filename = paste("fasta", ".fa", sep='')
	},
	content = function(file) {
	  selected <- get.contig.selected()
	  
	  if (any(selected$file.id%in%names(contig_data()))==FALSE) {
	    showModal(modalDialog(title = "ERROR!",
	                          paste0("Could not find a fasta file for: ",file),
	                          easyClose = T,
	                          footer = NULL))
	    return(NULL)
	  }
	  
	  apply(selected, 1, function(current){
	    
  	  contigid <- current['contig.id']
  	  fafile <- current['file.id']
  	  
  	  fasta <- contig_data()[[fafile]]$fasta
  	  idx <- contig_data()[[fafile]]$idx
  	  
  	  if (any(seqlevels(idx)==contigid)==FALSE) {
  	    showModal(modalDialog(title = "ERROR!",
  	                          paste0("Could not find contig: ",contigid," in your fasta file"),
  	                          easyClose = T,
  	                          footer = NULL))
  	    return(NULL)
  	  }
  	  
  	  loc <- idx[seqlevels(idx)==contigid]
  	  
  	  contig <- scanFa(fasta, param=loc)[[1]]
  	  contig <- DNAStringSet(contig)
  	  names(contig) <- contigid
  	  writeXStringSet(contig,file,append=T)
	  })
	}
)
