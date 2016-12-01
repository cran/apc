\name{apc-package}
\alias{apc-package}
\alias{apc}
\docType{package}
\title{Age-period-cohort analysis}
\description{The package includes functions for age-period-cohort analysis.  The statistical model is a generalized linear model (GLM)
allowing for age, period and cohort factors, or a sub-set of the factors.
The canonical parametrisation of Kuang, Nielsen and Nielsen (2008) is used. 
The outline of an analysis is described below.
}
\details{
	\tabular{ll}{
		Package: \tab apc\cr
		Type: \tab Package\cr
		Version: \tab 1.3\cr
		Date: \tab 2016-11-16\cr
		License: \tab GPL-3\cr
	}
	The apc package uses the canonical parameters suggested by
	Kuang, Nielsen and Nielsen (2008)
	and generalized by
	Nielsen (2014). 
	These evolve around the second differences of age, period and cohort factors as well as an three parameters (level and two slopes)
	for a linear plane.  The age, period and cohort factors themselves are not identifiable.  They could be ad hoc identified
	by associating the levels and two slopes to the age, period and cohort factors in a particular way.  This should be done
	with great care as such ad hoc identification easily masks which information is coming from the data and which information
	is coming from the choice of ad hoc identification scheme. An illustration is given below.
	A short description of the package can be found in
	Nielsen (2015).

	A formal analysis of the identification of the age-period-cohort model can be found in
	Nielsen and Nielsen (2014).
	Forecasting is not covered as by the package as yet, but discussion can be found in
	Kuang, Nielsen and Nielsen (2008b, 2011)
	and
	Martinez Miranda, Nielsen and Nielsen (2015).

	The apc package can be used as follows.
	\enumerate{
		\item
			Organize the data in as an \code{\link{apc.data.list}}.
			Data are included in matrix format.  Information needs to be given about the original data format.
			Optionally, information can be given about the labels for the time scales.
		\item
			Construct descriptive plots using \code{\link{apc.plot.data.all}}.
			This gives a series of descriptive plots.  The plots can be called individually through
			\enumerate{		
				\item
					Plot data sums using \code{\link{apc.plot.data.sums}}.
					Numerical values can be obtained through \code{\link{apc.data.sums}}.
				\item
				 	Sparsity plots of data using \code{\link{apc.plot.data.sparsity}}.
				\item	
				 	Plot data using all combinations of two time scales using \code{\link{apc.plot.data.within}}.
			}
		\item
			Get an deviance table for the age-period-cohort model through
			\code{\link{apc.fit.table}}.
		\item
			Estimate a particular (sub-model of) age-period-cohort model through
			\code{\link{apc.fit.model}}.
		\item
			Plot probability transforms of observed responses given fit using
			\code{\link{apc.plot.fit.pt}}.
		\item
			Plot estimated parameters through
			\code{\link{apc.plot.fit}}.
			Numerical values of certain transformations of the canonical parameter can be obtained through
			\code{\link{apc.identify}}.
		\item
			Recursive analysis can be done by selecting a subset of the observations through
			\code{\link{apc.data.list.subset}} and then repeating analysis.  This will reveal how sensitive
			the results are to particular age, period and cohort groups.
		\item
			Forecasting.  Some functions have been been added for forecasting in from a Poisson response-only model
			with an age-cohort parametrization
			\code{\link{apc.forecast.ac}}
			and with an age-period parametrization
			\code{\link{apc.forecast.ap}}.
			See also the overview on
			\code{\link{apc.forecast}}
	}
	Data examples include
	\enumerate{
		\item
			\code{\link{data.asbestos}}
			includes counts of deaths from mesothelioma in the UK.
			This dataset has no measure for exposure.
			It can be analysed using a Poisson model with an "APC" or an "AC" design.
			Source: Martinez Miranda, Nielsen and Nielsen (2015).
			Also used in Nielsen (2015).
		\item	
			\code{\link{data.Italian.bladder.cancer}}
			includes counts of deaths from bladder cancer in the Italy.
			This dataset includes a measure for exposure.
			It can be analysed using a Poisson model with an "APC" or an "AC" design.
			Source: Clayton and Schifflers (1987a).
		\item	
			\code{\link{data.Belgian.lung.cancer}}
			includes counts of deaths from lung cancer in the Belgium.
			This dataset includes a measure for exposure.
			It can be analysed using a Poisson model with an "APC", "AC", "AP" or "Ad" design.
			Source: Clayton and Schifflers (1987a).
		\item	
			\code{\link{data.Japanese.breast.cancer}}
			includes counts of deaths from breast cancer in the Japan.
			This dataset includes a measure for exposure.
			It can be analysed using a Poisson model with an "APC" design.
			Source: Clayton and Schifflers (1987b).
	}
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 29 Jan 2015 updated 7 Mar 2016.}
\references{
Clayton, D. and Schifflers, E. (1987a) Models for temperoral variation in cancer rates. I: age-period and age-cohort models. \emph{Statistics in Medicine} 6, 449-467.

Clayton, D. and Schifflers, E. (1987b) Models for temperoral variation in cancer rates. II: age-period-cohort models. \emph{Statistics in Medicine} 6, 469-481.

Kuang, D., Nielsen, B. and Nielsen, J.P. (2008a) Identification of the age-period-cohort model and the extended chain ladder model. Biometrika 95, 979-986. \emph{Download}: \href{http://biomet.oxfordjournals.org/cgi/reprint/95/4/979}{Article}; Earlier version \href{http://www.nuffield.ox.ac.uk/economics/papers/2007/w5/KuangNielsenNielsen07.pdf}{Nuffield DP}.

Kuang, D., Nielsen, B. and Nielsen, J.P. (2008b) Forecasting with the age-period-cohort model and the extended chain-ladder model. Biometrika 95, 987-991. \emph{Download}: \href{http://biomet.oxfordjournals.org/cgi/reprint/95/4/979}{Article}; Earlier version \href{http://www.nuffield.ox.ac.uk/economics/papers/2008/w9/KuangNielsenNielsen_Forecast.pdf}{Nuffield DP}.

Kuang, D., Nielsen, B. and Nielsen, J.P. (2011) Forecasting in an extended chain-ladder-type model. Journal of Risk and Insurance 78, 345-359. \emph{Download}: \href{http://dx.doi.org/10.1111/j.1539-6975.2010.01395.x}{Article}; Earlier version: \href{http://www.nuffield.ox.ac.uk/economics/papers/2010/w5/Forecast24jun10.pdf}{Nuffield DP}. 

Martinez Miranda, M.D., Nielsen, B. and Nielsen, J.P. (2015) Inference and forecasting in the age-period-cohort model with unknown exposure with an application to mesothelioma mortality. \emph{Journal of the Royal Statistical Society} A 178, 29-55. \emph{Download}: \href{http://www.nuffield.ox.ac.uk/economics/papers/2013/Asbestos8mar13.pdf}{Nuffield DP}.

Nielsen, B. (2015) apc: An R package for age-period-cohort analysis. \emph{R Journal} 7, 52-64. \emph{Download}: \href{https://journal.r-project.org/archive/2015-2/nielsen.pdf}{Open access}.

Nielsen, B. (2014) Deviance analysis of age-period-cohort models. \emph{Download}: \href{http://www.nuffield.ox.ac.uk/economics/papers/2014/apc_deviance.pdf}{Nuffield DP}.

Nielsen, B. and Nielsen, J.P. (2014) Identification and forecasting in mortality models. The Scientific World Journal. vol. 2014, Article ID 347043, 24 pages. \emph{Download}: \href{http://www.hindawi.com/journals/tswj/2014/347043}{Article}.
}
\seealso{
Vignettes are given on this
\href{http://users.ox.ac.uk/~nuff0078/apc/index.htm}{web page}.

Age-period-cohort analysis can alternatively be done by the package
\code{Epi}.
}
\examples{
########################
#	Belgian lung cancer

#######
#	1. Get apc.data.list
#	This is ready made.  For other data construct list using apc.data.list

data.list	<- data.Belgian.lung.cancer()
objects(data.list)
data.list

#######
#	2. Plot data
#	Plot all data.
#	Note a warning is produced because the defaults settings
#	lead to an unbalanced grouping of data.

apc.plot.data.all(data.list)

#	Or make individual plots.
#	Plot data sums.

graphics.off()
apc.plot.data.sums(data.list)

#	Plot sparsity to see where data are thin.
#	Plots are blank with default settings
#	... therefore change sparsity.limits.

graphics.off()
apc.plot.data.sparsity(data.list)
dev.new()
apc.plot.data.sparsity(data.list,sparsity.limits=c(5,10))

#	Plot data using different pairs of the three time scales.
#	This plot is done for mortality ratios.
#	All plots appear to have approximately parallel lines.
#	This indicates that interpretation should be done carefully.

graphics.off()
apc.plot.data.within.all.six(data.list,"m")

#######
#	3. Get a deviance table
#	Need to input distribution.
#   The table show that the sub-models "AC" and "Ad"
#	cannot be rejected relative to the unrestricted "APC" model

apc.fit.table(data.list,"poisson.dose.response")

#######
#	4. Estimate selected models
#	Consider "APC" and "Ad"
#	Consider also the sub-model "A", which is not supported by
#	the tests in the deviance table

fit.apc	<- apc.fit.model(data.list,"poisson.dose.response","APC")
fit.ad	<- apc.fit.model(data.list,"poisson.dose.response","Ad")
fit.a	<- apc.fit.model(data.list,"poisson.dose.response","A")

#	Get coefficients for canonical parameters through

fit.apc$coefficients.canonical
fit.ad$coefficients.canonical

#######
#	5. Residual analysis.
#	Plot estimators, probability transforms of responses given fit,
# 	residuals, fitted values, linear predictors, and data.
#	In probability transform plot:
#	Black circle are used for central part of distribution.
#	Triangles are used in tails, green/blue/red as responses are further in tail
#	No sign of mis-specification for "APC" and "Ad": there are many
#	black circles and only few coloured triangles.
#	In comparison the model "A" yields more extreme observations.
#	That model is not supported by the data.  
#	To get numerical values see apc.plot.fit.pt

graphics.off()
apc.plot.fit.all(fit.apc)
apc.plot.fit.all(fit.ad)
apc.plot.fit.all(fit.a)

#######
#	6. Plot estimated coefficients for sub models 
#	Consider "APC" and "Ad"
#	The first row of plots show double differences of paramters
#	The second row of plots shows level and slope determining a linear plane
#	The third row shows double sums of double differences,
#	all identified to be zero at the begining and at the end.
#	Thus the plots in third row must be interpreted jointly with those in the
#	second row.  The interpretation of the third row plots
#	is that they show deviations from linear trends.  The third row plots are
#	not invariant to changes to data array

graphics.off()
apc.plot.fit(fit.apc)
dev.new()
apc.plot.fit(fit.ad)
dev.new()
apc.plot.fit(fit.a)

#######
#	7. Recursive analysis
#	Cut the first period group and redo analysis

data.list.subset.1 <- apc.data.list.subset(data.list,0,0,1,0,0,0)
apc.fit.table(data.list.subset.1,"poisson.dose.response")

#######
#	8. Effect of ad hoc identification
#	At first a subset is chosen where youngest age and cohort groups
#	are truncated.  This way sparsity is eliminated
#	and ad hoc identification effects are dominated by estimation
#	uncertainty. Then consider
#	Plot 1: parameters estimated from data without first age groups
#	Plot 2: parameters estimated from all data
#	Note that estimates for double difference very similar.
#	Estimates for linear slopes are changed because the indices used
#	for parametrising these are changed
#	Estimates for detrended double sums of age and cohort double differences
#	are changed, because they rely on a particular ad hoc identifications
#	that have changed.  Nonetheless these plots are useful to evaulate
#	variation in time trends over and above linear trends.

graphics.off()
data.list	<- data.Belgian.lung.cancer()
data.list.subset <- apc.data.list.subset(data.list,2,0,0,0,0,0)
fit.apc		<- apc.fit.model(data.list,"poisson.dose.response","APC")
fit.apc.subset	<- apc.fit.model(data.list.subset,"poisson.dose.response","APC")
apc.plot.fit(fit.apc.subset,main.outer="1. Belgian lung cancer: cut first two age groups")
dev.new()
apc.plot.fit(fit.apc,main.outer="2. Belgian lung cancer data: all data")



}
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
\keyword{ hplot }
