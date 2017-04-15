#!/usr/bin/env Rscript

# options(warn=-1)

args = commandArgs(trailingOnly=TRUE) # TAKE ARGUMENTS FOLLOWING RSCRIPT CALL

directory = args[1] # directory containing the data
numberOfLoops = strtoi(args[2]) # number of trials to loop through

dataPaths = paste(directory, c('returns', 'assetConstraints', 'sectorConstraints', 'classConstraints'), sep='')

library(quadprog)

# FIXED OBJECT COMPUTATION

mean2 = function(x) {
mean(x, na.rm = TRUE)
}

returns = read.csv(paste(dataPaths[1], ".csv", sep=''), header=T) # import master dataframe from filepath first argument 

### COMPUTE FIXED PARAMETERS

assetNames = names(returns) # column headers; asset names
n = length(assetNames) # problem dimension; number of assets

X = scale(as.matrix(na.omit(returns)), center = TRUE, scale = FALSE) # centered feature-vector matrix
returns = apply(returns, 2, mean2) # overwrite master dataframe to save space
X = (X < 0) * X # downfall events
covar = t(X)%*%X # lower semi-covariance

EVD = eigen(covar) # spectral decomposition
minEig = min(EVD$values) # minimum eigenvalue
maxEig = max(EVD$values) # maximum eigenvalue

if (minEig < 0) { # if covariance is indefinite

	writeLines('\n\tSemi-Covariance is indefinite. Projecting into to semi-definite cone.')

	EVD$values = EVD$values - (EVD$values < 0)*(EVD$values - 1e-10) # spectral soft-thresholding
	covar = EVD$vectors%*%diag(EVD$values)%*%t(EVD$vectors) # reconstruct PSD covariance
	minEig = 1e-10 # new smallest eigenvalue

}

condNum = abs(maxEig/minEig) # condition number

if (condNum > 1e14) { # if the convariance is ill-conditioned

	writeLines('\tSemi-Covariance is ill-conditioned. Regularizing spectrum.')

	reg = maxEig/1e14 - minEig # regularizer enough for conditioning
	covar = covar + reg*diag(n) # then precondition
	remove(reg) # garbage collection

}

remove(X,EVD,minEig,maxEig,condNum) # garbage collection

### SIMPLEX COMPONENTS

d = matrix(rep(0,n),ncol=1) # auxiliary vector needed for quadprog to run
expansionPack = matrix(rep(1,n),nrow=1) # to copy b vec across rows

# CONSTRAINTS: SIMPLEX CONTRAINTS

Asimplex = rbind(
	returns, # must sum to 100%
	rep(1,n), # convex formulation constraint
	diag(n) # no shorts
	)

bsimplex = matrix(c(1,rep(0,n+1)),ncol=1) # right (in)equality vector

for (i in seq(1,numberOfLoops)) {

	args = paste(dataPaths[2:4], "_", toString(i), ".csv", sep='') # .../directory/data_iterationNumber.csv

	if (length(args)==0) {

			stop('ERROR: NO DATA FILEPATH PROVIDED.', call=TRUE)

		} else if (length(args)==3) {

			writeLines(paste('\n\tINITIALIZING ITERATION ',toString(i),'.',sep=''))

		} else {

			stop('ERROR: INCORRECT NUMBER OF FILEPATH ARGUMENTS PROVIDED.', call=TRUE)

		}

	### INITIALIZING QUADRATIC PROGRAMMING

	mEquality = 1 # stores number of equality constraints

	# CONSTRAINTS: ASSET ALLOCATION CONSTRAINTS

	assetConstraints = read.csv(args[1], header=T) # import asset constraints table
	attach(assetConstraints)
	# assetIndex = match(asset, assetNames)

	# minimum

	nonEmpty = asset[is.na(minimum) == FALSE] # select all assets with specified min
	m = length(nonEmpty) # number of min constraints

	if (m > 0) {

		b = minimum[is.na(minimum) == FALSE] # lower bounds

		AassetMin = -b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(1:m,match(nonEmpty,assetNames))) # relevant entries
		AassetMin[vMask] = AassetMin[vMask] + 1 # variable selection

	} else { AassetMin = NA; bassetMin = NA }

	# maximum

	nonEmpty = asset[is.na(maximum) == FALSE] # select all assets with specified max
	m = length(nonEmpty) # number of max constraints

	if (m > 0) {

		b = maximum[is.na(maximum) == FALSE] # upper bound

		AassetMax = b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(1:m,match(nonEmpty,assetNames))) # relevant entries
		AassetMax[vMask] = AassetMax[vMask] - 1 # variable selection

	} else { AassetMax = NA }

	# CONSTRAINTS: SECTOR LEVEL

	sectorConstraints = read.csv(args[2], header=T) # sector constraints table

	# # component-wise minimum

	# nonEmpty = sectorConstraints$sector[is.na(sectorConstraints$minimum) == FALSE] # sectors with minima
	# assetWithSector = subset(asset, sector%in%nonEmpty) # assets in sector with a min constraint
	# m = length(assetWithSector) # number of assets in indices with minima

	# AsectorMin = matrix(0, m, n) # for sector-level min
	# AsectorMin[cbind(1:m, match(assetWithSector, assetNames))] = -1 # select relevant assets

	# bsectorMin = -sectorConstraints$minimum[match(sector[match(assetWithSector, asset)],nonEmpty)]

	# # component-wise maximum

	# nonEmpty = sectorConstraints$sector[is.na(sectorConstraints$maximum) == FALSE] # sectors with maxima
	# assetWithSector = subset(asset, sector%in%nonEmpty) # assets in sector with a max constraint
	# m = length(assetWithSector) # number of assets in indices with maxima

	# AsectorMax = matrix(0, m, n) # for sector-level max
	# AsectorMax[cbind(1:m, match(assetWithSector, assetNames))] = 1 # select relevant assets

	# bsectorMax = sectorConstraints$maximum[match(sector[match(assetWithSector, asset)],nonEmpty)]

	# minimum

	nonEmpty = sectorConstraints$sector[is.na(sectorConstraints$minimum) == FALSE] # sectors with min specified
	assetWithSector = subset(asset, sector%in%nonEmpty) # assets in sector with a min constraint
	m = length(nonEmpty) # number of sectors with minima

	if (m > 0) {

		b = sectorConstraints$minimum[is.na(sectorConstraints$minimum) == FALSE] # lower bounds

		AsectorMin = -b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(match(sector[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AsectorMin[vMask] = AsectorMin[vMask] + 1 # variable selection

	} else { AsectorMin = NA }

	# maximum

	nonEmpty = sectorConstraints$sector[is.na(sectorConstraints$maximum) == FALSE] # sectors with max specified
	assetWithSector = subset(asset, sector%in%nonEmpty) # assets in sector with a max constraint
	m = length(nonEmpty) # number of sectors with max

	if (m > 0) {

		b = sectorConstraints$maximum[is.na(sectorConstraints$maximum) == FALSE] # upper bound

		AsectorMax = b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(match(sector[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AsectorMax[vMask] = AsectorMax[vMask] - 1 # select relevant assets

	} else { AsectorMax = NA }

	# sum

	nonEmpty = sectorConstraints$sector[is.na(sectorConstraints$summed) == FALSE] # sectors with sums specified
	assetWithSector = subset(asset, sector%in%nonEmpty) # assets in sector with a sum constraint
	m = length(nonEmpty) # number of sectors with sums

	if (m > 0) {

		b = sectorConstraints$summed[is.na(sectorConstraints$summed) == FALSE] # equality constraints

		AsectorSum = -b%*%expansionPack # for sector-level sum
		vMask = na.omit(cbind(match(sector[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AsectorSum[vMask] = AsectorSum[vMask] + 1 # select relevant assets

		mEquality = mEquality + length(b) # increase number of equality constraints
		
	} else { AsectorSum = NA }

	# CONSTRAINTS: CLASS LEVEL

	classConstraints = read.csv(args[3], header=T) # class constraints table

	# minimum

	nonEmpty = classConstraints$class[is.na(classConstraints$minimum) == FALSE] # classes with min specified
	assetWithClass = subset(asset, class%in%nonEmpty) # assets in class with a min constraint
	m = length(nonEmpty) # number of classes with minima

	if (m > 0) {

		b = classConstraints$minimum[is.na(classConstraints$minimum) == FALSE] # lower bound

		AclassMin = -b%*%expansionPack # initialize at zero # for class-level sum
		vMask = na.omit(cbind(match(class[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AclassMin[vMask] = AclassMin[vMask] + 1 # select relevant assets
		
	} else { AclassMin = NA }

	# maximum

	nonEmpty = classConstraints$class[is.na(classConstraints$maximum) == FALSE] # classes with max specified
	assetWithClass = subset(asset, class%in%nonEmpty) # assets in class with a max constraint
	m = length(nonEmpty) # number of classes with maxima

	if (m > 0) {

		b = classConstraints$maximum[is.na(classConstraints$maximum) == FALSE] # upper bound

		AclassMax = b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(match(class[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AclassMax[vMask] = AclassMax[vMask] - 1 # select relevant assets

	} else { AclassMax = NA }

	# sum

	nonEmpty = classConstraints$class[is.na(classConstraints$summed) == FALSE] # classes with sum specified
	assetWithClass = subset(asset, class%in%nonEmpty) # assets in class with a sum constraint
	m = length(nonEmpty) # number of classes with sums

	if (m > 0) {

		b = classConstraints$summed[is.na(classConstraints$summed) == FALSE] # equality constraint

		AclassSum = -b%*%expansionPack # initialize at zero
		vMask = na.omit(cbind(match(class[match(assetNames,asset)],nonEmpty), 1:n)) # relevant entries
		AclassSum[vMask] = AclassSum[vMask] + 1 # select relevant assets

		mEquality = mEquality + length(b) # increase number of equality constraints

	} else { AclassSum = NA }

	# CONSTRAINTS: CONCATENATION

	A = na.omit( # get rid of empty rows
		rbind( # row-wise concatenation of:
			AsectorSum,AclassSum, # equality constraints
			Asimplex, # simplex constraints
			AassetMin, AassetMax, # asset-level inequality constraints
			AsectorMin, AsectorMax, # sector-level inequality constraints
			AclassMin, AclassMax # class-level inequality constraints
		)
	)

	m = dim(A)[1] - (mEquality - 1) - length(bsimplex) # number of box inequlities

	b0 = c( # concatenation of:
		rep(0,mEquality-1), # equality bounds
		bsimplex, # simplex bounds
		rep(0,m) # asset-level, sector-level, and class-level bounds
	)

	# GARBAGE COLLECTION 

	detach(assetConstraints)

	remove(
		AassetMax, AassetMin, 
		AclassMax, AclassMin, AclassSum,
		args,
		AsectorMax, AsectorMin, AsectorSum,
		assetConstraints, assetWithClass, assetWithSector,
		b,
		classConstraints, sectorConstraints,
		m, nonEmpty, vMask
	)

	# COMPUTATION

	writeLines(paste('\tBEGINNING COMPUTATION ON ITERATION ',toString(i),'...',sep=''))

	qp = solve.QP(Dmat = covar, # covariance matrix in objective
		dvec = d, # equal to zero; must have so quadprog doesn't throw error
		Amat = t(A), bvec = b0, # constraints: non-negativity, weights sum to one
		meq = mEquality # only first m equations are equality (rest are inequality)
		)$solution # only store solution from the quadratic program

	writeLines('\tCOMPUTATION SUCCESSFUL!')

	### POST-PROCESSING + EXPORT

	weights = round(qp/sum(qp),7) # portfolio sum 500
	meanRet = sum(weights*returns) # expected returns
	volatility = sqrt(t(weights)%*%covar%*%weights) # portfolio covariance

	if (i == 1) { # initialize an empty df
		
		masterData = data.frame(
			fields=c('minAllocation','maxAllocation','expectedReturn','downsideRisk','sortinoRatio',assetNames),
			row.names='fields'
			)
	
	}

	masterData[paste('strategy',i, sep='_')] = c(min(weights),max(weights),round(meanRet,3),round(volatility,3),round(meanRet/volatility,3),weights)

}

write.csv(masterData, # export strategy comparison
file='./sortinoLoops.csv', # write it to current directory
row.names=TRUE) # include the row names

writeLines('\n\tSUCCESS. RESULTS WRITTEN TO CURRENT DIRECTORY AS sortinoLoops.csv!\n')
