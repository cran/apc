
\name{apc.identify}
\alias{apc.identify}
\title{Identification of double sums of canonical parameters}
\description{This computes double sums of canonical parameters along with detrended versions of these.  }
\usage{apc.identify(apc.fit.model)}
\arguments{
  \item{apc.fit.model}{List. See \code{\link{apc.fit.model}} for a description of the format.}
}
\details{Both types of double sums should be interpreted jointly with the linear trends.
The double sums have two zero values at the anchoring points defined in ??
The detrended sums start and end in zero.  They show the deviation of the time effects from linear trends.}
\value{
	\item{index.age.max		}{Vector.  Indices for age    parameters when using coefficients.ssdd or coefficients.detrend.  The length is two longer that that of \code{apc.model.fit$index.age} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.per.max		}{Vector.  Indices for period parameters when using coefficients.ssdd or coefficients.detrend.  The length is two longer that that of \code{apc.model.fit$index.per} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.coh.max		}{Vector.  Indices for cohort parameters when using coefficients.ssdd or coefficients.detrend.  The length is two longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{dates.max			}{Vector.  Indicates the dates for the parameters when using coefficients.ssdd or coefficients.detrend.  The length is six longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.}
	\item{index.age.sub		}{* Vector.  Indices for age    parameters when using coefficients.demean.  The length is two longer that that of \code{apc.model.fit$index.age} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.per.sub		}{* Vector.  Indices for period parameters when using coefficients.demean.  The length is two longer that that of \code{apc.model.fit$index.per} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.coh.sub		}{* Vector.  Indices for cohort parameters when using coefficients.demean.  The length is two longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{dates.sub			}{* Vector.  Indicates the dates for the parameters when using coefficients.demean.  The length is six longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.}
	\item{index.age.dif		}{* Vector.  Indices for age    parameters when using coefficients.dif.  The length is one longer that that of \code{apc.model.fit$index.age} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.per.dif		}{* Vector.  Indices for period parameters when using coefficients.dif.  The length is one  longer that that of \code{apc.model.fit$index.per} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{index.coh.dif		}{* Vector.  Indices for cohort parameters when using coefficients.dif.  The length is one  longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.  NULL if age    double differences are not estimated.}
	\item{dates.dif			}{* Vector.  Indicates the dates for the parameters when using coefficients.dif.  The length is three longer that that of \code{apc.model.fit$index.coh} if \code{model.design} is "APC.}
	\item{coefficients.ssdd		}{Matrix.  Coefficients of the double sum of double differences.  Normalised to be zero at two values chosen so age=cohort and period is at the minimal value.  For each parameter is reported coefficient, standard deviation, z-value, which is the ratio of those, and p-value.}
	\item{covariance.ssdd		}{Matrix.  Estimated covariance matrix for double sums.}
	\item{coefficients.detrend	}{Matrix.  Coefficients of the double sum of double differences.  Normalised to be zero for first and last value.  For each parameter is reported coefficient, standard deviation, z-value, which is the ratio of those, and p-value.}
	\item{covariance.detrend	}{Matrix.  Estimated covariance matrix for detrended double sums.}
	\item{coefficients.demean	}{* Matrix.  Coefficients of the sum of differences.  Normalised to be zero for first value.  Does not apply is design is "APC" For each parameter is reported coefficient, standard deviation, z-value, which is the ratio of those, and p-value.}
	\item{covariance.demean		}{* Matrix.  Estimated covariance matrix for demeaned sums.}
	\item{coefficients.dif		}{* Matrix.  Coefficients of the differences.  Does not apply is design is "APC" For each parameter is reported coefficient, standard deviation, z-value, which is the ratio of those, and p-value.}
	\item{covariance.dif		}{* Matrix.  Estimated covariance matrix for differences.}
}
\note{* indicates that values only implemented for designs "AC", "AP", "PC".

The differences are not identified for design "APC".  An arbitrary level can be moved between differences for
age, period and cohort.  

The differences are not identified for designs "Ad", "Pd", "Cd".  These models have two linear trends and one set of
double differences.  In the model "Ad", as an example, one linear trend will be associated with age, but it is
arbitrary whether the second linear trend should be associated with period or cohort.  The slope of the age trend
will depend on that arbitrary choice.  In turn the level of the age differences will be arbitrary.  
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 23 Aug 2014}
\examples{
########################
#	Belgian lung cancer
# 	first an example with APC design, note that demean and dif not defined.

data.list	<- data.Belgian.lung.cancer()

fit.apc	<- apc.fit.model(data.list,"poisson.dose.response","APC")
fit.apc$coefficients.canonical
id.apc	<- apc.identify(fit.apc)
id.apc$coefficients.ssdd
id.apc$coefficients.detrend
id.apc$coefficients.demean
id.apc$coefficients.dif

fit.ap	<- apc.fit.model(data.list,"poisson.dose.response","AP")
fit.ap$coefficients.canonical
id.ap	<- apc.identify(fit.ap)
id.ap$coefficients.ssdd
id.ap$coefficients.detrend
id.ap$coefficients.demean
id.ap$coefficients.dif

}