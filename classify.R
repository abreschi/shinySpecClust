suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(dtw))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(padr))
suppressPackageStartupMessages(library(imputeTS))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(pdist))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(xts))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mmand))

# ~~~~~~~
# OPTIONS
# ~~~~~~~

def_options = function() {
	option_list <- list(
	
		make_option(c("-i", "--input"), 
			help="File with CGM profile. 2 columns: 
			<date>, <value>. Has header [default=%default]"),
	
		make_option(c("-t", "--test_windows"), 
			help="File with windows. Provide in alternative to 
			--input. Format is windowId in col1 and raw values 
			in the other cols. Has NO header"),
	
		make_option(c("-P", "--parameters"), 
			help=".Rdata file with trained parameters."),
		
		make_option(c("-w", "--train_windows"), 
			help="File with training windows. 
			Has header, 1st col is window id, 
			remaining columns are CGM values"),
	
		make_option(c("-O", "--overlap"), default=0.25,
			help="Fraction of window overlap. For example 
			an overlap of 0.25 for 2.5 hour windows is
			37 mins, meaning a shift of 112 mins. [default=%default]"),
	
	#	make_option(c("-k", "--nb_clusters"), default=3,
	#		help="Number of desired clusters [default=%default]"),
	
		make_option(c("-o", "--output"), default="stdout",
			help="Output file name. [default=%default]"),
	
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
	return(opt)
}


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
get_na_stretches = function(x, minutes="90 min", 
		cgm_freq="5 min") {
	#minutes = 90
	#na_len = minutes/5
	na_len = duration(minutes)/duration(cgm_freq)
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
	smoothed = c(x[1:p], stats::filter(x, 
		smooth_weights/sum(smooth_weights), 
		sides=2)[(p+1):(n-p-1)], x[(n-p):n])
	return (smoothed)
}


smooth_WA = function(x) {	
	smooth_weights = c(1,2,4,8,16,24,16,8,4,2,1)
	smoothed = round(stats::filter(x, 
		smooth_weights/sum(smooth_weights), sides=2, 
		circular=T), 2)
	return (smoothed)
}


# preprocess
preprocess_cgm = function(m, gap='90 min', cgm_freq="5 min") {
	# m is a data.table with at least two columns
	# where the first is datetime string and 
	# the second is a numerical value

	# Convert to proper date time
	m[,1] = ymd_hms(m[[1]])
	
	# Thicken dates with 5 min frequency
	# Handle exception with already rounded datetimes
	thickened = tryCatch( { thicken( m, interval=cgm_freq)[,-1] }, 
		#error=function(c) m[,2:1] )
		error=function(c) {m[[1]] = m[[1]] + duration('1 sec');
		thicken( m, interval=cgm_freq)[,-1] } )
	thick_idx = ncol(thickened)

	# Make group for padding
	thickened$group = factor(c(0, 
		cumsum(diff(thickened[[thick_idx]]) > duration(gap))))

	# Pad dates
	padm = pad( thickened, group="group" )[,-c("group")]

	# Swap columns after padding
	#padm = padm[, c(2,1)]
	setcolorder(padm, colnames(padm)[c(
		thick_idx, 1:(thick_idx-1))])

	# Make sure glucose values are numeric
	padm[,2] = as.numeric(padm[[2]])

	## # Impute missing values up to <minutes>	
	## gluc_na_idx = get_na_stretches(padm[[2]], 
	## 	gap, cgm_freq=cgm_freq)

	# Impute missing values
	glucImp = na.interpolation(padm[[2]], option="stine")

	# Smooth with weigthed average
	smoothed = smooth_WA(glucImp)
	 
	## # Replace imputed values with NAs when 
	## # the stretch was too long
	## # Do this after smoothing
	## smoothed[gluc_na_idx] = NA

	# Assign smoothed values to padded datatable
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
	# Special case: all equal CGM values
	which_all_equal = which(apply(values, 1,
		function(x) sum(diff(x))) == 0)
	if (length(which_all_equal) > 0) {
		values[which_all_equal, 1] = values[which_all_equal, 1] + 1}
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
			sample(which(l==x), ceiling(mean(l==x)*N), 
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


prepare_test_set = function(cgm, param_list, window_overlap=0.25) {
	test = list()
	# Process query profile
	# ----------------------
	# Pad, fill and smooth cgm profile
	proc_cgm = preprocess_cgm(cgm)
	# Make windows of specified size and overlap
	# The size has to be fixed because it depends on the 
	# training set to compute the distance
	window_size = "150 mins"
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
	test_labels_na_omit = apply(tcrossprod(f, train$centers),
		1, which.max) 
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


prepare_test_windows = function(test_windows, param_list) {
	# Prepare test set from windows rather than from cgm profile
	test = list()
	test$windowID = test_windows[[1]]
	test$test = test_windows[,-1]
	# Scale baed on train parameters
	test$test = scale_windows(test$test, 
		param_list$train_mean, param_list$train_sd)
	# Smooth windows individually
	test$test = t(apply(test$test, 1, smooth_windows))
	# prepare labels for classification
	test$test_labels = rowSums(test$test)
	return(test)
}

classify_glucotype = function(cgm, train_windows, param_list, window_overlap) {
	# wrapper function for classifying windows
	test = prepare_test_set(cgm, param_list, window_overlap)
	train = prepare_training_set(test, train_windows, param_list)
	test = predict_windows(test, train)
	return (list('test'=test, 'train'=train))
}


classify_windows = function(cgm, train_windows, param_list, window_overlap) {
	# wrapper function for classifying windows
	res = classify_glucotype(cgm, train_windows,
		param_list, window_overlap)
	ts = reshape_test_windows(res$test, res$train)
	return (ts)
}

classify_windows2 = function(cgm, train_windows, param_list, window_overlap) {
	# wrapper function for predicting windows
	# same as classify_windows but the output has only the window labels
	res = classify_glucotype(cgm, train_windows,
		param_list, window_overlap)
	test = res$test
	train = res$train
	glucotypes = define_glucotypes(train)
	DT = data.table(test$cgm_w_windows)[windowPos==1,] %>% 
		.[,label:=glucotypes[test$test_labels[windowId]]] %>% 
		.[,c(4,6)]
	return (DT)
}

classify_test_windows = function(test_windows, train_windows, param_list) {
	# Classify windows without profile
	test = prepare_test_windows(test_windows, param_list)
	train = prepare_training_set(test, train_windows, param_list)
	test = predict_windows(test, train)
	DT = data.table(windowID = test$windowID, label=test$test_labels)
	return (DT)
}


cv = function(x) {
	sd(x,na.rm=T)/mean(x,n.rm=T)
}

rank_glucose_values = function(d) {
	# Finds most represented glucose value ranges,
	# out of 10 ranges, specific to the dataset.
	# Returns the ranking
	# --------------------------------------------
	# d is data.table with col1: timestamp, col2: glucose value
	# Already preprocessed
	cuts = as.numeric(cut(d[[2]], 10))
	gr = rank(table(cuts))[cuts]
	return(gr)
}

baseline_from_rank = function(d, gr) {
	# Get median of most frequently occurring glucose values
	# -------------------------------------------------------
	# d is data.table with col1: timestamp, col2: glucose value
	# Already preprocessed
	baseline = median(d[gr==9|gr==10][[2]])
	return(baseline)
}

get_baseline_rank = function(x) {
	# From unprocessed glucose values get baseline
	# based on most frequently occurring glucose values
	# -------------------------------------------------------
	# x is data.table with col1: timestamp, col2: glucose value	
	d = preprocess_cgm(x)
	gr = rank_glucose_values(d)
	baseline = baseline_from_rank(d, gr)
	return (baseline)
}

smooth_cq = function(cq) {
	cqs = as.integer(rollapply(zoo(cq), 5, 
		median, na.rm=T, na.pad=T, partial=T,
		coredata = TRUE))
	close_gaps = closing(is.na(cq), c(1,1,1,1,1))
	cqs = replace(cqs, 
		which(as.logical(close_gaps)), NA)
	return(cqs)
}

decreasing_cq = function(windows) {
	# Returns negative trend of windows
	#gradient_qs = apply(windows, 1, 
	#	function(x) mean(diff(x)))
	gradient_qs = apply(windows, 1, function(x) 
		quantile(sign(diff(x)), 0.75, na.rm=T))
	# 1 when trend is negative, 0 otherwise
	decr_bin = replace(rep(0, nrow(windows)), 
		gradient_qs<=0, 1)
	return (decr_bin)
}


cq_gap_wilcox = function(windows, x, y) {
	pvalue = mapply(function(x,y) {
		wilcox.test(unlist(windows[x,]), unlist(windows[y,]), 
			'greater')$p.value
		}, x, y)
	return(pvalue)
}

fill_gaps = function(cqs, windows) {
# -- Fill the gaps --
	# Find stretches of consecutive windows in the same 30% quantile
	stretch = rle(cqs==1 & !is.na(cqs))
	# Total length of adjacent low variance windows without gap
	ll = rollsum(stretch$lengths[stretch$values], 2)
	# Total length of adjacent low variance windows including gap
	llgap = rollsum(stretch$lengths, 3)[stretch$values]
	# Gap lengths
	gaps = stretch$lengths[which(stretch$values)[1:(sum(stretch$values)-1)] +1]
	# Indeces of gaps
	gaps_idx = which(stretch$values)[1:(sum(stretch$values)-1)] +1
	
	# Corresponding window indeces for gaps, before, and after gaps	
	gaps_wins = mapply(seq, c(0,cumsum(stretch$lengths))[gaps_idx]+1,  
		cumsum(stretch$lengths)[gaps_idx])
	gaps_before_wins = mapply(seq, c(0,cumsum(stretch$lengths))[gaps_idx-1]+1,  
		cumsum(stretch$lengths)[gaps_idx-1])
	gaps_after_wins = mapply(seq, c(0,cumsum(stretch$lengths))[gaps_idx+1]+1,  
		cumsum(stretch$lengths)[gaps_idx+1])
	
	# Test before and after 
	before_greater = cq_gap_wilcox(windows, gaps_before_wins, gaps_wins)
	after_greater = cq_gap_wilcox(windows, gaps_after_wins, gaps_wins)
	
	# Get the index of gap to fill
	gaps_fill_idx = which(gaps <= 5 & ll/llgap <= 2/3 & before_greater <= 0.05 & after_greater <= 0.05)
	
	# Fill gaps in smoothed quantiles
	cqsf = replace(cqs, unlist(gaps_wins[gaps_fill_idx]), 1)
	
	return (cqsf)
}

longest_stretch = function(x, fill=3) {
	# Initialize resulting vector
	xo = rep(0, length(x))
	## Close gaps if values are spearated by <fill>
	#x = closing(x, rep(1, fill))
	# Count consecutive stretches
	wr = rle(as.numeric(x == 1 & !is.na(x)))
	# Find relative id of longest stretch of flat windows
	ll = which.max(wr$lengths[wr$values == 1])
	# Make sure there is at least one flat period
	if ( length(ll) == 0 ) {return(xo)}
	# Find absolute id of longest stretch of flat windows
	lla = which(wr$values == 1)[ll]
	# Check if the longest stretch is at the beginning or end of window
	if ( lla == 1 | lla == length(wr$length) ) {return(xo)}
	# Check if longest stretch is in the middle of the window but the window starts and/or end with a stretch
	if ( (ll > 1 & ll < sum(wr$values)) & (wr$values[1] == 1 | wr$values[length(wr$values)] == 1) ) {return(xo)}
	#cumsum(wr$lengths)[wr$values == 0][ll] > cumsum(wr$lengths)[wr$values == 1][ll]
	# Assign 1s to resulting vector only for longest stretch
	xo[ max(1, cumsum(wr$lengths)[lla-1]) : cumsum(wr$lengths)[lla] ] <- 1
#	print(cumsum(wr$lengths)[max(1, lla-1)] : cumsum(wr$lengths)[lla] )
	return(xo)
}

sleep_periods = function(d) {
	# Find periods of stable glucose values from 
	# processed cgm data
	# ------------------------------------------
	# Extract column names
	timecol = colnames(d)[1]
	valuecol = colnames(d)[2]
	# Generate the windows
	windows = make_windows(d, '2.5 hours', 0.75)
	wins = data.table(annotate_cgm_windows(d, "2.5 hours", 0.75))
	# Assign windows to quantiles based on their coefficient of variation
	win_cv = apply(windows, 1, cv)
	cq = as.numeric(cut(win_cv, c(0, quantile(win_cv, 
		probs=c(3:10)/10, na.rm=T))))
	## Add decreasing windows to low variability windows, even 
	## if they have high variability
	#decr_bin = decreasing_cq(windows)
	#cq[decr_bin == 1 & !is.na(cq)] = 1
	wins$cq = cq[wins$windowId]
	# Smooth quantiles - TODO: improve padding
	cqs = smooth_cq(cq) 
	wins$cqs = cqs[wins[["windowId"]]]
	# Fill the gaps
	cqsf = fill_gaps(cqs, windows)
	# Reassign quantiles after gap filling
	wins$cqsf = cqsf[wins[["windowId"]]]
	# Find the sleep periods by shifting 24 hour-window
	ann_wins = shift_24_hours(wins, d)
	
	gp = ggplot(data.frame(cv=win_cv)) + 
		geom_histogram(aes(cv), fill='salmon') + 
		geom_vline(data=data.frame(cv=c(0, 
		quantile(win_cv, probs=1:10/10, na.rm=T))), 
		aes(xintercept=cv)) + 
		geom_vline(xintercept=quantile(win_cv, probs=0.3, na.rm=T), size=3)

	gp = ggplot(
			cbind( 
				wins, 
				days = floor_date(wins[[timecol]], 'days'),
				xend = wins[[timecol]] + hours(2) + minutes(30)
			),
			aes_string(timecol, valuecol)) + 
		geom_line() + 
		geom_point(aes(color=as.factor(cqs))) + 
		scale_color_manual(values=rainbow(10)) + 
		facet_wrap(~days, scales='free_x') + 
		geom_vline(data=wins[wins[, un:=length(unique(cqs))==1, 
			by=c(timecol)][["un"]]][cqs==1,], 
			aes_string(xintercept=timecol), alpha=0.2, size=2) +
		geom_vline(data=data.table(wins)[cqs==1], aes_string(xintercept=timecol), 
			alpha=0.2, size=2, color='pink') 
	#	geom_hline(yintercept=92) +
	#	geom_hline(yintercept=100)
	#ggsave('tmp.pdf', plot=gp, h=15, w=20)
	# Find the longest periods of low variation per day
	return(ann_wins)
}


shift_24_hours = function(wins, d) {
	# Move a 24 hour-windows every 6 hours to find
	# the longest stretches of stable glucose every 
	# daily cycle.
	# ---------------------------------------------
	# Extract column names
	timecol = colnames(d)[1]
	valuecol = colnames(d)[2]
	# Parameters for windows
	win_size = 24
	win_shift = 6
	# Select only one cqsf for each timepoint (the minimum) 
	ann_wins = wins[, list(cqsf=as.numeric(min(cqsf)==1)), 
		by=c(timecol, valuecol)]
	ann_wins[, win_idx:= 1+ as.numeric(seconds(interval(first(get(timecol)), 
		get(timecol)))) %/% (win_size*3600)]
	ann_wins[, shift_idx:= 1+ as.numeric(seconds(interval(first(get(timecol)), 
		get(timecol)))) %/% (win_shift*3600)]
	ann_wins$flats = rep(0, nrow(ann_wins))
	# Find the longest stretch within each window
	for ( i in 0:(win_size/win_shift -1) ) {
		ann_wins$factor = ann_wins$win_idx + 
			as.numeric((ann_wins$shift_idx-1) %% (win_size/win_shift) >= i) -1; 
		flats_i = ann_wins[, longest_stretch(cqsf), by="factor"][["V1"]]; 
		ann_wins$flats = as.numeric( ann_wins$flats | flats_i) 
		# Plot glucose values with shaded shifted windows of 24 hours
		gp = ggplot(cbind(ann_wins, col=ann_wins[["factor"]]%%2==1, 
				days=floor_date(ann_wins[[colnames(d)[1]]], "days")), 
			aes_string(colnames(d)[1], colnames(d)[2])) + 
			facet_wrap(~days, scales="free_x") + 
			geom_vline(aes_string(xintercept=timecol, alpha="flats"), color='pink', size=2) + 
			geom_line() + 
			geom_vline(aes_string(xintercept=timecol, color="col"), alpha=0.1, size=2) + 
			scale_color_manual(values=c("orange", "blue"))
	}
	return(ann_wins)
} 

ann_wins_to_intervals = function(ann_wins, d) {
	# Extract intervals of sleep or stable glucose
	# from annotated timepoints
	# --------------------------------------------
	# Extract column names
	timecol = colnames(d)[1]
	valuecol = colnames(d)[2]
	# Find contiguous timepoints of stable glucose
	# or estimated sleep
	stretches = rle(
		replace(
			replace(
				replace(ann_wins$cqsf, is.na(ann_wins$cqsf), 0), 
				ann_wins$cqsf == 1, "stable"), 
			ann_wins$flats == 1, "sleep"
		)
	)
	indices = cumsum(stretches$lengths)
	intervals_indices = which(stretches$values != "0")
	intervals_ends = indices[intervals_indices]
	intervals_starts = indices[intervals_indices - 1] +1
	# Special case the interval is at the beginning of recordings
	if (intervals_indices[1] == 1 ) {
		intervals_starts = c(1, intervals_starts) }
	intervals = data.table(
		start = ann_wins[intervals_starts, ][[timecol]], 
		end = ann_wins[intervals_ends, ][[timecol]],
		type = stretches$values[intervals_indices],
		median = mapply( function(x,y) {
			median(ann_wins[x:y,][[valuecol]])}, 
				intervals_starts, intervals_ends
		)
	)
	intervals = adjust_sleep_intervals(intervals)
	return(intervals)
}



adjust_sleep_intervals = function(intervals) {
	# Consider sleep only if it overlaps
	# most recurrent skeep hour bin
	hours_sleep = unlist(apply(
		intervals[type=="sleep"], 1, function(x) 
		hour(seq(ymd_hms(x[1]), ymd_hms(x[2]), 
		by='hours'))))
	most_freq_hour = as.numeric(names(
		sort(table(hours_sleep), dec=T)[1]))
	a = ymd_hms(intervals[[2]])
	b = ymd_hms(intervals[[2]])
	hour(a) = most_freq_hour
	hour(b) = most_freq_hour
	minute(a) = 0
	minute(b) = 0
	most_freq_int = interval(a-hours(1),b+hours(1))
	int = interval(intervals[[1]], intervals[[2]])
	intervals[["adjusted"]] = "stable"
	adjusted_sleep = which(mapply(function(x,y) 
		as.logical(lubridate::intersect(x,y)), 
		most_freq_int, interval(intervals[[1]], 
		intervals[[2]])))
	intervals[["adjusted"]][adjusted_sleep] = "sleep"
	return (intervals)
}

get_baselines = function(d, ann_wins) {
	# Compute baselines obtained with 
	# different methods
	# ------------------------------------
	# Get baseline based on ranking
	gr = rank_glucose_values(d)
	median_rank = baseline_from_rank(d, gr) 
	
	# Medians - baseline
	median_grand = median(d[[2]], na.rm=T)
	
	# Baseline from sleep periods
	median_cq = median(ann_wins[flats==1,][[colnames(d)[2]]])	
	
	# Compile the data.table with baselines
	medians = data.table(
		method = c(
			'median_rank',
			'median_grand',
			'median_cq'
		),
		value = c(
			median_rank,
			median_grand,
			median_cq
		)
	)
	return(medians)
}	


# ~~~~~~~~~~~
# BEGIN
# ~~~~~~~~~~~

args = commandArgs(FALSE)
fileArg = args[grep("--file", args)]
script = ""
if ( length(fileArg) != 0 ){
	script = strsplit(strsplit(fileArg, 
		"=")[[1]][[2]], "/")[[1]]
	script = script[length(script)]
}


if(length(args)!=0 & script == "classify.R") {
	
	opt = def_options()

	print("Running")

	if (!is.null(opt$input) & !is.null(opt$test_windows)) {
		cat("Error: Provide cgm profile OR windows\n")
		q(save='no')
	}

	# Input files
	cgmF = opt$input; if(!is.null(opt$input)) {if( opt$input == "stdin") {
		cgmF = 'file:///dev/stdin'}}
	testF = opt$test_windows; if(!is.null(opt$test_windows)) { 
		if (opt$test_windows == "stdin") {
		testF = 'file:///dev/stdin'}}
	windowsF = opt$train_windows
	paramF = opt$parameters
	
	# Output
	outF = ifelse(opt$output == "stdout", 
		"", opt$output)
	# Window overlap
	window_overlap = opt$overlap
	
#	# Input variables for debugging
#	paramF = "train.overlap_37+window_2.5.params.Rdata"
#	cgmF = "cgm.test.tsv"
#	windowsF = "0_preprocessing/raw/rawDexcomSeries+overlap_37+window_2.5+user_all"
	
	# read input files
	if (!is.null(cgmF)) {cgm = fread(cgmF)}
	if (!is.null(testF)) {test_windows = fread(testF, h=F)}
	train_windows = fread(windowsF)
	load(paramF)
	
	# Classify windows from profile
	if(!is.null(cgmF)) {
		pred = classify_windows2(cgm, train_windows, 
		param_list, window_overlap)	
	}
	
	# Classify windows directly (not profile)
	if (!is.null(testF)) {
		pred = classify_test_windows(test_windows,
		train_windows, param_list)
	}


	# Write predictions to file
	write.table(pred, outF, sep='\t', quote=F, row.names=F)

}

