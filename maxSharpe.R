#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # TAKE ARGUMENTS FOLLOWING RSCRIPT CALL

if (length(args)==0) {

		stop('ERROR: NO DATA FILEPATH PROVIDED.', call=TRUE)

	} else if (length(args)==1) {

		writeLines('\n\tINITIALIZING MAX SHARPE RATIO QUADRATIC PROGRAM.\n')

	} else {

		stop('ERROR: CAN ONLY SUPPLY ONE FILEPATH ARGUMENT.', call=TRUE)

	}


### INITIALIZATION

library(quadprog)

mean2 = function(x) {
	mean(x, na.rm = TRUE)
}

returns = read.csv(args[1], header=T) # import master dataframe from filepath first argument 

### COMPUTE FIXED PARAMETERS

assetNames = names(returns) # column headers; asset names
n = length(assetNames) # problem dimension; number of assets

covar = cov(returns, use="complete.obs", method="pearson") # robust covariance; casewise deletion of missing data
eig = eigen(covar)$values # eigenvalues
minEig = min(eig) # smallest eigenvalue
maxEig = max(eig) # largest eigenvalue
condNum = abs(maxEig/minEig) # condition number

if (minEig < 0) { # if covariance is indefinite

	covar = covar - 2*minEig*diag(n) # push to semidefinite cone
	minEig = abs(minEig) # minimum eigenvalue is now positive

}

if (condNum > 1e13) { # if the convariance is ill-conditioned

	reg = maxEig/1e13 - minEig
	covar = covar + reg*diag(n) # then precondition

}

returns = apply(returns, 2, mean2) # overwrite master dataframe to save space

remove(eig,minEig,maxEig,condNum,reg) # garbage collection

### QUADRATIC PROGRAMMING

A = rbind(returns, rep(1,n), diag(n)) # left constraints equation matrix
b0 = matrix(c(1,rep(0,n+1)),ncol=1) # right (in)equality vector
d = matrix(rep(0,n),ncol=1) # auxiliary vector needed for quadprog to run

writeLines('\tBEGINNING COMPUTATION...\n')

qp = solve.QP(Dmat = covar, # covariance matrix in objective
	dvec = d, # equal to zero; must have so quadprog doesn't throw error
	Amat = t(A), bvec = b0, # constraints: non-negativity, weights sum to one
	meq = 1 # only first equation is equality (rest are inequality)
	)$solution # only store solution from the quadratic program

### POST-PROCESSING + EXPORT

weights = round(qp/sum(qp),7) # portfolio sum 500

exportData = data.frame(asset=assetNames, weight=weights) # table with allocations

write.csv(exportData, # export the weights as csv
	file='./maxSharpe.csv', # write it to current directory
	row.names=FALSE) # do not include the row index

meanRet = sum(weights*returns) # expected returns
volatility = sqrt(t(weights)%*%covar%*%weights) # portfolio covariance
sharpeRatio = meanRet/volatility

write.csv(data.frame(meanRet, volatility, sharpeRatio), # dataframe with portfolio stats
	file='./portfolioStats.csv', # write it to current directory
	row.names=FALSE) # do not include the row index

writeLines('\tSUCCESS. RESULTS WRITTEN TO CURRENT DIRECTORY AS maxSharpe.csv AND portfolioStats.csv!\n')
