library(shiny)
library(DT)
library(dygraphs)
library(ggplot2)
library(lubridate)
library(xts)
library(scales)


# Import functions for window classification
source("classify.R")

# Import spectral files and parameters
paramF = "train.overlap_37+window_2.5.params.Rdata"
load(paramF)
windowsF= "rawDexcomSeries+overlap_37+window_2.5+user_all"
train_windows = fread(windowsF)

# Other functions for data formatting
df_to_xts = function(df) {
	ts = xts(df[,-1], order.by=ymd_hms(df[,1]))
	return(ts)
}

# Dygraph parameters
alpha = 0.5
colors = paste("rgba(", apply( rbind( 
	col2rgb(rev(hue_pal()(3))), 
	alpha=rep(alpha, 3)), 2, paste,
	collapse=","), ")"
)
axisTitle = "Glucose (mg/dL)"
# Glucotypes
glucotypes = c("low", "moderate", "severe")

# Spinner image
spinnerUI = tags$img(src="loading_spinner.gif", 
	alt="spinner", 
	id="spinner",
	style ="position:absolute; margin: auto; z-index:2"
)


#add_spinner = "console.log('hello')"
#
#add_spinner = "$(document).on('shiny:inputchanged', function(event) {
#				if (event.name === 'glucotype') {
#					var elem = document.createElement('img');
#					elem.setAttribute('src', 'loading_spinner.gif');
#					elem.setAttribute('alt', 'spinner');
#					elem.setAttribute('id', 'spinner');
#					elem.style.position = 'absolute';
#					elem.style.margin = '0 auto';
#					elem.style['z-index'] = '2';
#					document.getElementById('dygraphContainer').prepend(elem);
#					console.log(document.getElementById('dygraphContainer'));
#				}
#				
#			})"




# ~~~~~~~~~
# UI
# ~~~~~~~~~

ui <- fluidPage(

#	tags$head(HTML("<link rel='stylesheet' src='normalize.css'>")),
	tags$head(tags$style(HTML("
		* {
			font-size: 20px;
		}
		p {
			text-align: justify;
		}
		h1 {
			font-size: 4rem;
		}
		body {
			margin-left: 2rem;
			margin-right: 4rem;
			max-width: 1800px;
		}
		.action-button {
			font-size: 1em;
		}
		.input-group * {
			font-size: 1em;
			height: auto;
		}
	"))),

	title = "CGM viewer",
	HTML("
		<div style='display: flex'> 
			<h1 style='align-self: center; flex-grow: 1'>CGM viewer</h1>
			<img style='width:auto; max-height: 80px; margin: 3rem' 
				src='Stanford_Medicine_logo-web-CS.trim.png' alt='logo'>
		</div>
	"),

	
	# Controls
	fluidRow(

		column(4, 
			HTML("<p>
			Visualize your CGM profile and discover your glucotype! 
			This webpage provides an interface to the glucotype classification
			described in <a href='#'>Hall et al</a>. 
			To find out about your glucotype please upload a tab-delimited
			file with two columns, containing the time of measurement and
			the glucose concentration, following the format in 
			<a href='sample.cgm.tsv' target='_blank'>this sample file</a>.
			</p> 
			")
		),
		
		column(3,
			style="margin-left: 4rem;",
			fileInput(inputId = "cgmF",
				"Upload CGM profile",
				multiple = TRUE,
				accept = c("text/csv",
					"text/comma-separated-values,text/plain",
					".csv", ".tsv", ".txt"
				)
			),

			tags$label("or    "),
			
			actionButton(inputId = "showSample",
				"Visualize sample profile"
			)
		),

		column(2, 

			actionButton(inputId = "glucotype",
				"Classify glucotype",
			),
			actionButton(inputId = "reset",
				"Clear"
			),
			br(),
			uiOutput(outputId = "verdict")
		)	
	),
	
	fluidRow(
		column(12,
			tags$div(id = "dygraphContainer",
				style = "position: relative; text-align: center; margin: 1rem auto;",
				tags$div(
					dygraphOutput(outputId = "cgmProfile"),
					br(),
					uiOutput(outputId = "dygraphDesc")
				)
			)
		)
	),

	br(),

	fluidRow(
		column(5, 
			div(DT::dataTableOutput("cgmTable"),
				style = "overflow-y:scroll; max-height:300px; max-width: 100%")
		),
		column(1),
		column(3,
			tableOutput(outputId = "glucotypeTable")
		)
	)

#	tags$head(HTML(sprintf("<script>%s</script>", add_spinner)))
)



# ~~~~~~~~~~~~~~~~
# SERVER
# ~~~~~~~~~~~~~~~~
# Define server logic
server <- function(input, output, session) {

	print("Restart")

	# Read 	
	read_cgmF = function(f) {
		df = read.csv(f,
			header = TRUE,
			sep = "\t",
			quote = ""
		)
		df[[1]] = ymd_hms(df[[1]])
		return (df)
	}

	v = reactiveValues(data=NULL);

	observeEvent(input$cgmF, {
		req(input[['cgmF']]);
		v$data = read_cgmF(input$cgmF$datapath)
		print(head(v$data))
	})

	observeEvent(input[["showSample"]], {
		v$data = read_cgmF("www/sample.cgm.tsv")
		print(head(v$data))
	})


	# table of raw cgm values
	output[['cgmTable']] = DT::renderDataTable({
		DT::datatable(v$data, 
			options = list())
	})

	# dygraph of raw cgm values
	dygraph_cgm = renderDygraph({
		if (is.null(v$data)) return()
		df = v$data
		dygraph(df_to_xts(df), ylab=axisTitle) %>%
		dyOptions(axisLabelWidth=90)
	})
	output[["cgmProfile"]] = dygraph_cgm
	# Description for dygraph navigation commands
	output[["dygraphDesc"]] = renderUI({
		if (is.null(v$data)) return()
		tags$p("Zoom: click-drag; Pan: shift-click-drag; 
		Restore zoom level: double-click", 
		style="text-align: center; color: grey")
	})



	# Calculate glucotype
	observeEvent(input$glucotype, {
#		cgm = data();
		cgm = v$data
		print("Compute glucotype")
		print(head(cgm))
		insertUI("#dygraphContainer", where="afterBegin",
			tags$div("Computing.... Please be patient", id="spinnerBg", 
				style="position: absolute; vertical-align: middle;
				background-color:rgba(200,200,200,0.4); 
				z-index:10; width:100%; height:100%"), 
			immediate=TRUE, session=session);
#		insertUI("#dygraphContainer", where="afterBegin", 
#			tags$div(id="spinnerBg", style="position:absolute;
#				background-color:rgba(200,200,200,0.4); 
#				z-index:10; width:100%; height:100%"), 
#			immediate=TRUE, session=session);
#		insertUI("#dygraphContainer", where="afterBegin", 
#			spinnerUI, immediate=TRUE, session=session);

		cachedF = "glucotypes.df.tsv"
		if (!file.exists(cachedF)) {
			#print(head(preprocess_cgm(cgm)));
			df = classify_windows(cgm, train_windows, param_list, 0.25)
#			Have to think of a better way to cache dataframe
#			write.zoo(df, file = cachedF, quote=F, sep='\t')
		} else {
			df = fread(cachedF)
			df = as.xts.data.table(df[,Index:=ymd_hms(Index)])
		}
		# remove loading spinner at the end of computation
#		observe({removeUI("#spinner"); })
		observe({removeUI("#spinnerBg"); })
		# Plot colored dygraph
		output[["cgmProfile"]] = renderDygraph({
			plot_colors = colors[match(colnames(df), glucotypes)]
			print(head(df))
			dygraph(df) %>%
			dyAxis("y", label = axisTitle) %>%
			dyOptions(axisLabelWidth=90) %>%
			dyOptions(colors=plot_colors, drawPoints=T, pointSize=2)
		})
		# Table with glucotype frequencies
		freq = data.frame("Fraction of time" = 
				colMeans(!is.na(df)), check.names=F)
		output[["glucotypeTable"]] = renderTable(
			rownames = TRUE, {
			freq
		})
		# Print glucotype on top of the page
		output[["verdict"]] = renderUI({
			gt = rownames(freq)[which.max(freq[[1]])]
			gt_color = rev(hue_pal()(3))[match(gt, glucotypes)]
			h1(gt, style=sprintf("color: %s", gt_color))
		})
	})


	# Remove glucotypes and table
	observeEvent(input$reset, {
		output[["cgmProfile"]] = dygraph_cgm
		output[["glucotypeTable"]] = renderTable({})
	})

}

shinyOptions = list(
	port = 8000, 
	launch.browser = T
)

shinyApp(ui = ui, server = server, options = shinyOptions)
