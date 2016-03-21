\name{apc.forecast.ap}
\alias{apc.forecast.ap}
\title{Forecast for Poisson response model with AP structure.}
\description{Computes forecasts for a model with AP structure.  Forecasts of the linear predictor are given for all models.
For the forecast the period parameters need to be extrapolated.
The extrapolation method has to chosen so as not to introduce an identification problem, see
Kuang, Nielsen and Nielsen (2008b,2011).
Two extrapolation methods are possible: "I0" and "I1".
Those paper 
}
\usage{apc.forecast.ap(apc.fit,extrapolation.type="I0",suppress.warning=TRUE)}
\arguments{
  \item{apc.fit}{List. Output from \code{\link{apc.fit.model}}.
  					Note: apc.fit.model should be run with AP structure so that 
					\code{apc.fit$model.design=="AP"}.
					Only point forecasts of the linear predictor are provided.}
  \item{extrapolation.type}{Character.  Choices for extrapolating the differenced period parameter ("Delta.beta_per").
  							Default is "I0".
  							\describe{
  								\item{"I0"}{extrapolates the first out-of-sample differenced period parameter by the average of
								cumulated sums of the in-sample estimated differenced period parameters.
								The subsequent out-of-sample differenced period parameters are zero.}
								\item{"I1"}{extrapolates all out-of-sample differenced period parameters by zero.}
							}
							Both methods are invariant to ad hoc identification of the implied period time effect, by
							following the ideas put forward in
							Kuang, Nielsen and Nielsen (2008b).
							Internally, the extrapolation is done as follows.
							The estimated differenced period parameters are found from
							"apc.fit$coefficients.canonical" using
							\code{\link{apc.identify}} with \code{type="dif"}.
							These imply period time effects by ad hoc identification:
							choose an arbitrary value for the first period time effect
							and add partial sums of the differenced period parameter. Fit a time series model:
							an intercept model with "I0" and a random walk model for "I1".
							Then extrapolate and take differences.  These extrapolation methods are invariant to the
							actual choice of the arbitrary value for the first period time effect.}
	\item{suppress.warning}{Logical. If true, suppresses warnings from \code{\link{apc.data.list.subset}}, which is called internally. Default is "TRUE".}							
}
\value{
  \item{trap.linear.predictors.forecast}{Matrix. Includes estimates and point forecasts of linear predictor. Forecasts in lower right triangle. Trapezoid format.}
  \item{index.trap.J}{Matrix. age-coh coordinates for forecast area. Similar structure to
  					\code{index.trap} in \code{apc.index}, see \code{\link{apc.get.index}}.}
}
\references{
Kuang, D., Nielsen, B. and Nielsen, J.P. (2008b) Forecasting with the age-period-cohort model and the extended chain-ladder model. Biometrika 95, 987-991. \emph{Download}: \href{http://biomet.oxfordjournals.org/cgi/reprint/95/4/987}{Article}; Earlier version \href{http://www.nuffield.ox.ac.uk/economics/papers/2008/w9/KuangNielsenNielsen_Forecast.pdf}{Nuffield DP}.

Kuang, D., Nielsen B. and Nielsen J.P. (2011) Forecasting in an extended chain-ladder-type model. Journal of Risk and Insurance 78, 345-359. \emph{Download}: \href{http://dx.doi.org/10.1111/j.1539-6975.2010.01395.x}{Article}; Earlier version: \href{http://www.nuffield.ox.ac.uk/economics/papers/2010/w5/Forecast24jun10.pdf}{Nuffield DP}. 
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 2 Mar 2016}