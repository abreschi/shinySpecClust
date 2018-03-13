#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtw))
suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(optparse))



# ~~~~~~~
# OPTIONS
# ~~~~~~~

option_list <- list(

	make_option(c("-i", "--input"), default="stdin",
		help="File with preprocessed windows [default=%default]"),
	
	make_option(c("-k", "--nb_clusters"), default=3,
		help="Number of desired clusters [default=%default]"),

	make_option(c("-d", "--distances"), 
		help="File with pre-computed distance matrix"),

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

# Smoothing function, to perform after scaling
smooth_windows = function(x) {
	smooth_weights = c(1,2,4,8,16,8,4,2,1)
	p = (length(smooth_weights)-1)/2
	n = length(x)
	smoothed = c(x[1:p], filter(x, 
		smooth_weights/sum(smooth_weights), 
		sides=2)[(p+1):(n-p-1)], x[(n-p):n])
	return (smoothed)
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

get_labels = function(Y) {
	eigDiscrete = .discretisation(Y)
	eigDiscrete = eigDiscrete$discrete
	labels = apply(eigDiscrete, 1, which.max)
	centers = as.matrix(aggregate(Y, mean, by=list(labels))[,-1])
	kM = list(cluster=labels, centers=centers)
	return (kM)
}

# Get training parameters. These should not be
# recomputed every time to make predictions
train_parameters = function(train, nb_clusters, distances=NULL) {
	# train is a datatable with n timeseries as rows
	# each column represents a different time point
	
	# List of parameters to save in a separate parameter file
	param_list = list()	

	# Get grand mean and sd
	train_mean = mean(as.matrix(train))
	train_sd = sd(as.matrix(train))
	param_list[['train_mean']] = train_mean
	param_list[['train_sd']] = train_sd
	
	# Scale training matrix
	train_z = (train - train_mean)/train_sd
	# Smooth
	train_sm = data.table(t(apply(train_z, 1, smooth_windows)))
	
	# Compute normalization factors for DTW
	if (opt$verbose) {cat("Computing CID-DTW normalization factors...")}
	train_CE = apply(train_sm, 1, function(x) sqrt(sum(diff(x)^2)))
	train_CF = outer(train_CE, train_CE, Vectorize( function(x,y) 
		return(max(x,y)/min(x,y))))
	param_list[['train_CE']] = train_CE
	if (opt$verbose) {cat("DONE\n")}
	# Read precomputed distances
	if (!is.null(distances)) {
		cat("Skipping distance computation\n")
		train_dist = distances
	} else {
		# Compute dtw distance for training
		if (opt$verbose) {cat("Computing training DTW distances...")}
		shakoechiba_window_size = round(ncol(train)*0.1, 0)
		train_dist = as.matrix(dist(train_sm, method="DTW", window.type = "sakoechiba", 
			window.size = shakoechiba_window_size, step.pattern = symmetric2)) * train_CF
		rm(train_CF)
		if (opt$verbose) {cat("DONE\n")}
	}
	
	# Number of neighbors for affinity matrix
	K = round(sqrt(nrow(train)))
	# Compute sigma for training
	train_dist_sorted = as.matrix(t(apply(train_dist, 2, sort)))
	train_means = rowMeans(train_dist_sorted[, 1:K + 1]) + .Machine$double.eps
	rm(train_dist_sorted)
	param_list[["train_means"]] = train_means

	# compute affinity matrix for training
	if (opt$verbose) {cat("Computing training affinity matrix...")}
	A = affMatrix( train_dist, train_means, 0.5)
	if (opt$verbose) {cat("DONE\n")}

	# Compute diagonal matrix for spectral clustering
#	diag(A) = 0
	if (opt$verbose) {cat("Computing unnormalized Laplacian...")}
	D = diag(rowSums(A))
	L = D - A
	rm(A)
	if (opt$verbose) {cat("DONE\n")}

	if (opt$verbose) {cat("Computing normalized Laplacian...")}
	L = diag(1/sqrt(diag(D))) %*% L %*% diag(1/sqrt(diag(D)))
	rm(D)
	if (opt$verbose) {cat("DONE\n")}

	if (opt$verbose) {cat("Computing eigen values/vectors...")}
	eig = eigen(L, symmetric=T)
	rm(L)
	res = sort(abs(eig$values), index.return = TRUE)
	X = eig$vectors[, res$ix[1:nb_clusters]]
	eig_values = eig$values
	rm(eig)
	param_list[['eig_values']] = eig_values
	if (opt$verbose) {cat("DONE\n")}

	if (opt$verbose) {cat("Normalizing eigen vectors...")}
	# Normalized eigen vectors 
	# keep only first <nb_clusters> eigen vectors
	Y = X/sqrt(rowSums(X^2))
	param_list[["Y"]] = Y
	if (opt$verbose) {cat("DONE\n")}

	if (opt$verbose) {cat("Computing centroids...")}
	# Compute centroids based on pre-identified labels
	kM = get_labels(Y)
#	centers = aggregate(Y, by=list(labels[[2]]), FUN=mean)
#	centers = as.matrix(centers[,-1])
	param_list[['centers']] = kM$centers
	param_list[['labels']] = kM$cluster
	if (opt$verbose) {cat("DONE\n")}

	return(param_list)
}


# ~~~~~~
# BEGIN
# ~~~~~~

# Read options
#labelsF = opt$labels
f = opt$input
nb_clusters = opt$nb_clusters
outputF = opt$output
distancesF = opt$distances
distances = NULL


## debugging variables
#f = "0_preporocessing/raw/rawDexcomSeries+overlap_75+window_2.5+user_all"
#nb_clusters = 3
#labelsF = "train.labels.sample.tsv"
#outputF = "train.param_list.Rdata"

# Open files and store tables
#labels = fread(labelsF, sep='\t') 
train = fread(f, sep='\t')
windows = train[[1]]
train = train[,-1]
# Read distances if provided
if (!is.null(distancesF)) {
	distances = fread(input = sprintf('zcat < %s', distancesF), h=T)
	distances = as.matrix(as.dist(distances[,-1]))
}

# Get training parameters
param_list = train_parameters( train, nb_clusters, distances)
# Save them to file
save(param_list, list="param_list", file=outputF)

q(save='no')
