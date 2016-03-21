\name{apc.forecast}
\alias{apc.forecast}
\title{Forecasts from age-period-cohort models.}
\description{In general forecasts from age-period-cohort models require extrapolation of the estimated parameters.
	This has to be done without introducing identifications problems, see
	Kuang, Nielsen and Nielsen (2008b,2011).
	There are many different possibilities for extrapolation for the different sub-models.
	The extrapolation results in point forecasts.
	Distribution forecasts should be build on top of these, see
	Martinez Miranda, Nielsen and Nielsen (2015).
	At present two experimental functions
	\code{\link{apc.forecast.ac}}
	and
	\code{\link{apc.forecast.ap}}
	are available.
}
\references{
Kuang, D., Nielsen, B. and Nielsen, J.P. (2008b) Forecasting with the age-period-cohort model and the extended chain-ladder model. Biometrika 95, 987-991. \emph{Download}: \href{http://biomet.oxfordjournals.org/cgi/reprint/95/4/987}{Article}; Earlier version \href{http://www.nuffield.ox.ac.uk/economics/papers/2008/w9/KuangNielsenNielsen_Forecast.pdf}{Nuffield DP}.

Kuang, D., Nielsen B. and Nielsen J.P. (2011) Forecasting in an extended chain-ladder-type model. Journal of Risk and Insurance 78, 345-359. \emph{Download}: \href{http://dx.doi.org/10.1111/j.1539-6975.2010.01395.x}{Article}; Earlier version: \href{http://www.nuffield.ox.ac.uk/economics/papers/2010/w5/Forecast24jun10.pdf}{Nuffield DP}. 

Martinez Miranda, M.D., Nielsen, B. and Nielsen, J.P. (2015) Inference and forecasting in the age-period-cohort model with unknown exposure with an application to mesothelioma mortality. \emph{Journal of the Royal Statistical Society} A 178, 29-55. \emph{Download}: \href{http://www.nuffield.ox.ac.uk/economics/papers/2013/Asbestos8mar13.pdf}{Nuffield DP}.
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 1 Feb 2016}