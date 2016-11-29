\name{covEB}
\alias{covEB}
\title{
Empirical Bayes estimate of block diagonal correlation matrix
}
\description{
Assuming a block diagonal structure of the correlation matrix, the function calculates and empirical Bayes estimate of the original covariance matrix. The aim is to reduce false discovery rates by pooling information on the levels of correlations between elements in the same blocks. The blocks can have different levels of true correlation between them. The algorithm searches the space of possible correlation values and estimates the final correlation by and average of all non zero estimates. 
}
\usage{
covEB(Covmat, delta = 0.1, shift = 0.05, cutoff = NULL, startlambda = 0.5,n)
}
\arguments{
  \item{Covmat}{
A sample covariance matrix. Must be positive semi-definite, is converted to a correlation matrix in the algorithm.
}
  \item{delta}{
This parameter gives the size of the range of correlation values that are considered in one iteration. The larger the value the more the values will be pooled to the overall correlation of the full matrix. Default 0.1
}
  \item{shift}{
The shift parameter determines how much the lower bound of the range of correlation values is moved by through each iteration. The larger the value the fewer iterations. Default 0.05.
}
  \item{cutoff}{
Optional - If the average correlations for a given block are below this value then the prior for this block is set to be the independence prior.
}

  \item{startlambda}{
This determines a lower noise level for the correlations, for example, it defaults to 0.5. This means that all correlation values below 0.5 will be set to zero (assumed to be noise).
}

\item{n}{
The number of samples (replicates) used to calculate Covmat}
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

	test<-covEB(samplecov,delta=0.05,shift=0.025,startlambda=0.4,n=500)
}


\keyword{ correlation }
\keyword{ empirical }
\keyword{ bayes }
