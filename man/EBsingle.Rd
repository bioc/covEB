\name{EBsingle}
\alias{EBsingle}
\title{
Empirical Bayes estimate of block diagonal correlation matrix
}
\description{
Assuming a block diagonal structure of the correlation matrix, the function calculates and empirical Bayes estimate of the original covariance matrix. The algorithm assumes a single block diagonal prior with one shrinkage threshold determining significance of correlations. 
}
\usage{
EBsingle(Covmat, startlambda = 0.5, n, happrox=FALSE,groups=NULL,aic=TRUE,data=NULL)
}
\arguments{
  \item{Covmat}{
A sample covariance matrix. Must be positive semi-definite, is converted to a correlation matrix in the algorithm.
}


  \item{startlambda}{
This determines a lower noise level for the correlations, for example, it defaults to 0.5. This means that all correlation values below 0.5 will be set to zero (assumed to be noise). This will be used if groups parameter is not specified.
}

\item{n}{
The number of samples (replicates) used to calculate Covmat}

\item{happrox}{
Logical indicating whether or not to use a hypergeometric distribution estimation of the correlations. The alternative is to take the sample average of the correlations. FALSE by default.
}
\item{groups}{
	Optional list of groups of variables from Covmat that are used to define the block diagonal structure of the prior.
}
\item{aic}{
	Optional logical parameter indicating whether or not the AIC metric should be used to select the threshold parameter from the data. This will replace the value of startlambda. To use the AIC option, the data must also be provided
}
\item{data}{
	Optional data matrix, required when AIC is to be used. The matrix should have observations in rows and variables in the columns.
}
}

\value{
Returns the estimated correlation matrix.
}
\references{
Champion, C. J. (2003). Empirical Bayesian estimation of normal variances and covariances. Journal of Multivariate Analysis, 87(1), 60-79
}
\author{
C. Pacini
}

\examples{


	
	sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)

	samplecov<-cov(x)

	test<-EBsingle(samplecov,startlambda=0.4,n=500)
}


\keyword{ correlation }
\keyword{ empirical }
\keyword{ bayes }
