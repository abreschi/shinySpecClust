library(data.table)
library(lubridate)
library(dtw)
library(ggplot2)
library(proxy)
library(padr)
library(imputeTS)
library(zoo)
library(pdist)
library(optparse)
library(xts)

# ~~~~~~~
# OPTIONS
# ~~~~~~~

option_list <- list(

	make_option(c("-i", "--input"), default="stdin",
		help="File with CGM profile. 2 columns: 
		<date>, <value>. Has header [default=%default]"),

	make_option(c("-P", "--parameters"), 
		help=".Rdata file with trained parameters."),
	
	make_option(c("-w", "--train_windows"), 
		help="File with training windows. 
		Has header, 1st col is window id, 
		remaining columns are CGM values"),

	make_option(c("-k", "--nb_clusters"), default=3,
		help="Number of desired clusters [default=%default]"),

	make_option(c("-l", "--labels"),
		help="Two columns file with id and label. Rows
		should be in the same order as input"),

	make_option(c("-o", "--output"), default="train.param_list.Rdata",
		help="Output file name. Should be .Rdata [default=%default]"),

	make_option(c("-v", "--verbose"), default=FALSE, action="store_true",
		help="Verbose output")

)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = ""
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


# ~~~~~~~~~
# FUNCTIONS
# ~~~~~~~~~

# Split string with time interval into a list
make_interval_list_from_string <- function(interval_string) {
  interval_split <- strsplit(interval_string, " ")[[1]]
  if (length(interval_split) == 1) {
    return(list(interval = interval_split,
                step     = 1))
  } else {
    return(list(interval = interval_split[2],
                step     = as.numeric(interval_split[1])))
  }
}

# Interpret time units from different formats
uniform_interval_name <- function(interval) {
  if (interval %in% c("y", "ye", "yea", "years")) {
    interval <- "years"
  } else if (interval %in% c("q", "qu", "qua", "quar", "quart", "quarte", "quarters")){
    interval <- "quarters"
  } else if (interval %in% c("m", "mo", "mon", "mont", "months")) {
    interval <- "months"
  } else if (interval %in% c("w", "we", "wee", "weeks")){
    interval <- "weeks"
  } else if (interval %in% c("d", "da", "days")) {
    interval <- "days"
  } else if (interval %in% c("h", "ho", "hou", "hours")) {
    interval <- "hours"
  } else if (interval %in% c("mi", "mins")) {
    interval <- "mins"
  } else if (interval %in% c("s", "se", "secs")) {
    interval <- "secs"
  }
  return(interval)
}

# Convert interval to seconds
interval_to_seconds = function(interval) {
	interval_split = make_interval_list_from_string(interval)
	interval_split$interval = uniform_interval_name(interval_split$interval)
	tdiff = as.difftime(interval_split$step, 
		units=interval_split$interval)
	secs = as.numeric(tdiff, units='secs')
	return (secs)
}

# Convert interval in seconds to 5 minutes span
seconds_to_5_mins = function(secs) {
	mins = round(secs/60 / 5)
	return (mins)
}

# Convert interval to 5 minutes span
interval_to_5_mins = function(interval) {
	secs = interval_to_seconds(interval)
	mins = seconds_to_5_mins(secs)
	return (mins)
}

# create the distance metric
cid_dtw_dist <- function(time.seriesQ, time.seriesC, shakoechiba.window.size=2){
          CE_Q = sqrt(sum(diff(time.seriesQ)^2))
          CE_C = sqrt(sum(diff(time.seriesC)^2))
          CF = max(CE_Q,CE_C)/min(CE_Q,CE_C)
          dic_dtw = CF*dtw(time.seriesQ, time.seriesC, 
          window.type = "sakoechiba", window.size = shakoechiba.window.size,
          step.pattern = symmetric2,
          distance.only = T)$distance
      return(dic_dtw)
}

# add new distance metric to distance db (proxy package)
add_cid_dtw_to_db = function() {
	pr_DB$set_entry(FUN = cid_dtw_dist, names = c('cid-dtw'))
}

# 
get_na_stretches = function(x, minutes) {
	minutes = 90
	na_len = minutes/5
	gluc_na = rle(is.na(x))
	gluc_na$values = gluc_na$values & gluc_na$lengths >= na_len
	gluc_na_idx = which(inverse.rle(gluc_na))
	# Return which indexes have a stretch of NAs
	# longer than <minutes>
	return(gluc_na_idx)
}

smooth_windows = function(x) {
	smooth_weights = c(1,2,4,8,16,8,4,2,1)
	p = (length(smooth_weights)-1)/2
	n = length(x)
	smoothed = c(x[1:p], filter(x, 
		smooth_weights/sum(smooth_weights), 
		sides=2)[(p+1):(n-p-1)], x[(n-p):n])
	return (smoothed)
}


smooth_WA = function(x) {	
	smooth_weights = c(1,2,4,8,16,24,16,8,4,2,1)
	smoothed = round(filter(x, 
		smooth_weights/sum(smooth_weights), sides=2, 
		circular=T), 2)
	return (smoothed)
}


# preprocess
preprocess_cgm = function(m) {
	# Convert to proper date time
	m[,1] = ymd_hms(m[[1]])
	
	# Pad dates
	padm = pad( thicken( m, interval="5 min" )[,-1])
	# Swap columns after padding
	padm = padm[, c(2,1)]

	# Impute missing values up to <minutes>	
	gluc_na_idx = get_na_stretches(padm[[2]], minutes)

	# Impute missing values
	glucImp = na.interpolation(padm[[2]], option="stine")

	# Smooth with weigthed average
	smoothed = smooth_WA(glucImp)
	 
	# Replace imputed values with NAs when 
	# the stretch was too long
	# Do this after smoothing
	smoothed[gluc_na_idx] = NA

	# Assign smoothed values to padded datable
	padm[[2]] = smoothed

	return (padm)
}

# Determine shift based on percentage overlap
window_overlap_to_shift = function(window_size, w_overlap) {
	# overlap is given as fraction
	window_size = interval_to_seconds(window_size)
	# Shift in seconds
	shift = window_size - window_size * w_overlap
	shift = seconds_to_5_mins(shift)
	return(shift)
}

# Make overlapping windows of CGM data
make_windows = function(m, window_size, w_overlap) {
	size = interval_to_5_mins(window_size)
	shift = window_overlap_to_shift(window_size, w_overlap)
	values = as.data.table(rollapply(zoo(m[[2]]), 
		size, function(x) return(x), by=shift))
	window_starts = as.numeric(rollapply(zoo(m[[1]]), 
		size, function(x) return(x[1]), by=shift))
#	print(as_datetime(window_starts))
	return(values)
}

# Annotate the cgm profile with the windows indexes
annotate_cgm_windows = function(m, window_size, w_overlap) {
	size = interval_to_5_mins(window_size)
	shift = window_overlap_to_shift(window_size, w_overlap)
	window_idx = melt(
		rollapply(1:nrow(m), size, 
			function(x) return(x), by=shift), 
		varnames=c("windowId", "windowPos"), 
		value.name="cgmIdx"
	)
	# Merge window indices with cgm profile
	df = merge(window_idx, m, 
		by.x="cgmIdx", by.y="row.names")
	return(df)
}


# Scale windows based on pre-computed 
# mean and standard deviation
scale_windows = function(m, mean, sd) {
	return ( (m-mean)/sd )
}



# Compute affinity matrix from distance
affMatrix = function(train_dist, train_means, alpha) {
	# train_dist is a distance matrix
	# is symmetric in the case of the training
	# for testing the size mxn, where m is 
	# the number of testing examples and n is 
	# the number of examples in the training
	# Have to decide if keeping columns as training
	# or reshaping the matrix
	return (t(   # Have to transpose for sapply
		sapply(1:nrow(train_dist), function(x) {
			Sig = (train_means[x]+train_means)/3 + 
				train_dist[x,]/3 + .Machine$double.eps; 
			return( dnorm(train_dist[x,], 0, alpha*Sig))
			}
		)
	)	
)}



# Get representative points from kmeans
sample_training = function(Y, kM, N) {
	# Y are the standardized eigenVectors
	# kM are the results from k-means
	# N is the number of final examples
	# ------------------------------
	# Compute the distance between 
	# eigenVectors and each center
	set.seed(1234)
	dk = as.matrix(pdist(Y, kM$centers))
	l = kM$cluster
	idx = unlist( 
		sapply(unique(l), function(x) 
			sample(which(l==x), round(mean(l==x)*N), 
			prob=dk[cbind(1:length(l), l)][l==x]) 
		)
	)
	return (idx)
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}


# Apply .discretisation function from SNFtool package
get_labels = function(Y) {
	eigDiscrete = .discretisation(Y)
	eigDiscrete = eigDiscrete$discrete
	labels = apply(eigDiscrete, 1, which.max)
	centers = as.matrix(aggregate(Y, mean, by=list(labels))[,-1])
	kM = list(cluster=labels, centers=centers)
	return (kM)
}


prepare_test_set = function(cgm, param_list) {
	test = list()
	# Process query profile
	# ----------------------
	# Pad, fill and smooth cgm profile
	proc_cgm = preprocess_cgm(cgm)
	# Make windows of specified size and overlap
	# The size has to be fixed because it depends on the 
	# training set to compute the distance
	window_size = "150 mins"
	window_overlap = 0.25
	test$test = make_windows(proc_cgm, window_size, window_overlap)
	# Annotate cgm profile with window id
	test$cgm_w_windows = annotate_cgm_windows(proc_cgm, window_size, window_overlap)
	# Make a vector with NAs
	test$test_labels = rowSums(test$test)
	# Remove windows with NAs
	test$test = na.omit(test$test)
	# Scale windows
	test$test = scale_windows(test$test, 
		param_list$train_mean, param_list$train_sd)
	return(test)
}


prepare_training_set = function(test, train_windows, param_list) {
	train_param = list()
	# Get training parameters
	# ------------------------
	Y = param_list$Y
	# Get cluster labels
	kM = get_labels(Y)
	#set.seed(123); kM = kmeans(Y, centers=3, iter.max=20)
	# Choose a smaller sample of representative points
	nb_examples = max(200, nrow(test$test)) 
	train_idx = sample_training(Y, kM, nb_examples)
	train_param$train_Y = Y[train_idx,]
	# Get training windows
	train_param$train = scale_windows(train_windows[train_idx,-1], 
		param_list$train_mean, param_list$train_sd)
	train_param$train_means = param_list$train_means[train_idx]
	train_param$train_CE = param_list$train_CE[train_idx]
	train_param$train_labels = kM$cluster[train_idx]
	train_param$centers = kM$centers
	return(train_param)
}


predict_windows = function(test, train) { 
	# Classification of windows
	# --------------------------
	# Compute normalization factors for DTW
	test_CE = apply(test$test, 1, function(x) sqrt(sum(diff(x)^2)))
	# --- Check how long this will take for larger training ---
	test_CF = outer(test_CE, train$train_CE, Vectorize( function(x,y) 
		return(max(x,y)/min(x,y))))
	# This is very slow, 
	# consider moving to python fastDTW
	shakoechiba_window_size = round(ncol(train$train)*0.1, 0)
	test_dist = dtwDist(test$test, train$train, window.type = "sakoechiba", 
		window.size = shakoechiba_window_size, step.pattern = symmetric2) * test_CF
	
	# Compute affinity matrix for testing
	B = affMatrix(test_dist, train$train_means, 0.5)
	
	# Project on eigenfunctions
	f = B %*% train$train_Y
	# Normalize eigenfunctions
	f = f/sqrt(rowSums(f^2))
	
	# This computes the distance to the k-means 
	# centers on the transformed space
	test_labels_na_omit = apply(f, 1, function(x) { 
		which.max(tcrossprod(x, train$centers)) } )
	# Restore NAs
	test$test_labels[!is.na(test$test_labels)] = test_labels_na_omit 
	#df = merge(melt(as.matrix(test)), as.data.frame(test_labels), by.x="Var1", by.y='row.names')
	#ggplot(na.omit(df), aes(Var2, value)) + geom_line(aes(group=Var1, color=test_labels)) + facet_wrap(~test_labels)
	return(test)
}



define_glucotypes = function(train) {
	# Give labels to clusters based on the mean 
	# of the training windows
	levels = c("low", "moderate", "severe")
	means = apply(train$train, 1, mean)
	glucotypes = factor(levels[order(
		aggregate(means, median, 
		by=list(train$train_labels))$x)], 
		levels = levels
	)
	return(glucotypes)
}


reshape_test_windows = function(test, train) {
	glucotypes = define_glucotypes(train)
	df = cbind(test$cgm_w_windows[,c(4,5)],
		glucotype = glucotypes[ test$test_labels[ 
			test$cgm_w_windows[["windowId"]] ]
		]
	)
	cast_formula = as.formula( paste( 
		colnames(df)[1], "~glucotype"))
	ts = as.xts.data.table(dcast.data.table( 
		data.table(na.omit(df)), 
		cast_formula, 
		value.var="GlucoseValue", mean)
	)	
	return(ts)
}


classify_windows = function(cgm, train_windows, param_list) {
	# wrapper function for classifying windows
	test = prepare_test_set(cgm, param_list)
	train = prepare_training_set(test, train_windows, param_list)
	test = predict_windows(test, train)
	ts = reshape_test_windows(test, train)
}



# ~~~~~~~~~~~
# BEGIN
# ~~~~~~~~~~~

args = commandArgs(TRUE)

if(length(args)!=0) {

	print("Running")

	# Input files
	cgmF = opt$input
	windowsF = opt$train_windows
	paramF = opt$parameters
	
	# Input variables for debugging
	paramF = "train.overlap_37+window_2.5.params.Rdata"
	cgmF = "cgm.test.tsv"
	windowsF = "0_preprocessing/raw/rawDexcomSeries+overlap_37+window_2.5+user_all"
	
	# read input files
	cgm = fread(cgmF)
	train_windows = fread(windowsF)
	load(paramF)
	
	# Classify windows
	ts = classify_windows(cgm, train_windows, param_list)

}


