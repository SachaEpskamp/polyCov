\name{polyCov}
\alias{polyCov}

\title{
Polychoric Covariances
}
\description{
This function can be used to estimate a polychoric (or tetrachoric)  correlation 
matrix assuming standardized latent items. If thresholds are given, it estimates 
polychoric covariances, means and variances.
}
\usage{
polyCov(x, y = NULL, thresholds = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix or dataframe containing ordinal data, or a vector containing ordinal data.
}
  \item{y}{
A vector containing ordinal data. Only if x is also a vector.
}
  \item{thresholds}{
An optional list containing thresholds for each item.
}
  \item{\dots}{
Arguments sent to \code{optimize}
}
}
\details{
If 'thresholds' is not assigned, then the latent means and variances of each item 
are fixed to 0 and 1 respectively resulting in a correlation matrix. If `thresholds'
is assigned then the means and variances are estimated for each item and a covariance 
matrix is returned. For dichotomous items the variance is fixed to 1.

The polychoric correlations are estimated in a two-step procedure (Luo, 2011). first 
the thresholds are obtained from the marginal distribution of each variable and then 
the covariances are optimized to the frequency table.
}
\value{
	\item{covmat}{The estimated covariance matrix.}
	\item{thresholds}{A list containing the estimated thresholds.}
	\item{means}{A vector of estimated means.}
	\item{vars}{A vector of estimated variances.}
}
\references{
Luo, H. (2011) Some Aspects on Confirmatory Factor Analysis of Ordinal Variables 
and Generating Non-normal Data. Retrieved from http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-149423
}
\author{
Sacha Epskamp
}
\examples{
# Simulate data:
dat <- rnorm(1000) + matrix(rnorm(1000*5),,5)

# Real correlation matrix:
cor(dat)

# Transform data in ordinal scale:
datOrd <- ifelse(dat < -1,1,ifelse(dat>1,3,2))

# Polychoric correlations:
polyCov(datOrd)

}
\keyword{Ordinal data}
\keyword{Covariances}
\keyword{Correlations}
