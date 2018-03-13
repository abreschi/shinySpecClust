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


	# App title ----
	titlePanel("CGM viewer"),
	
	br(),

	# Controls
	fluidRow(

		column(4,
			fileInput(inputId = "cgmF", 
				"Upload CGM profile",
				multiple = TRUE,
				accept = c("text/csv",
					"text/comma-separated-values,text/plain",
					".csv", ".tsv", ".txt"
				)
			)
		),

		column(4,
			actionButton(inputId = "glucotype",
				"Classify glucotype"
			),
			actionButton(inputId = "reset",
				"Clear"
			)
		)	
	),
	
	fluidRow(
		column(10,
			tags$div(id = "dygraphContainer",
				style = "position: relative; text-align: center",
				tags$div(
					dygraphOutput(outputId = "cgmProfile")
				)
			)
		),
		column(2,
			plotOutput("plot")
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
	

	read_cgmF = reactive ({
		req(input[['cgmF']]);
		df = read.csv(input$cgmF$datapath,
			header = TRUE,
			sep = "\t",
			quote = ""
		)
		df[[1]] = ymd_hms(df[[1]])
		return(df)
	})

	# table of raw cgm values
	output[['cgmTable']] = DT::renderDataTable({
		DT::datatable(read_cgmF(), 
			options = list())
	})

	# dygraph of raw cgm values
	dygraph_cgm = renderDygraph({
		df = read_cgmF()
		print(head(df))
		dygraph(df_to_xts(df)) %>%
		dyAxis("y", label = axisTitle)
	})
	output[["cgmProfile"]] = dygraph_cgm;

	observeEvent(input$cgmF, {
		output[["cgmProfile"]] = dygraph_cgm
		output[["glucotypeTable"]] = renderTable({})
	})


	# Calculate glucotype
	observeEvent(input$glucotype, {
		cgm = read_cgmF();
		insertUI("#dygraphContainer", where="afterBegin", 
			tags$div(id="spinnerBg", style="position:absolute;
				background-color:rgba(200,200,200,0.4); 
				z-index:10; width:100%; height:100%"), 
			immediate=TRUE, session=session);
		insertUI("#dygraphContainer", where="afterBegin", 
			spinnerUI, immediate=TRUE, session=session);

		cachedF = "glucotypes.df.tsv"
		if (!file.exists(cachedF)) {
			#print(head(preprocess_cgm(cgm)));
			df = classify_windows(cgm, train_windows, param_list)
#			Have to think of a better way to cache dataframe
#			write.zoo(df, file = cachedF, quote=F, sep='\t')
		} else {
			df = fread(cachedF)
			df = as.xts.data.table(df[,Index:=ymd_hms(Index)])
		}
		# remove loading spinner at the end of computation
		observe({removeUI("#spinner"); })
		observe({removeUI("#spinnerBg"); })
		# Plot colored dygraph
		output[["cgmProfile"]] = renderDygraph({
			plot_colors = colors[match(colnames(df), glucotypes)]
			print(head(df))
			dygraph(df) %>%
			dyAxis("y", label = axisTitle) %>%
			dyOptions(colors=plot_colors, drawPoints=T, pointSize=2)
		})
		# Table with glucotype frequencies
		output[["glucotypeTable"]] = renderTable(
			rownames = TRUE, {
			data.frame("Fraction of time" = 
				colMeans(!is.na(df)), check.names=F)
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
