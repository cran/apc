\name{apc.plot.fit}
\alias{apc.plot.fit}
\title{Plots of apc estimates}
\description{Functions to plot the apc estimates found by \code{\link{apc.fit.model}}.  The function apc.plot.fit detects the type of
\code{model.design} and \code{model.family} from the fit values and makes appropriate plots.

Depending on the \code{model.design} the plot has up to 9 sub plots.
The type of these can be chosen using \code{type}

\bold{Model designs of any type.}
If \code{type} is "detrend" or "sum.sum"
the canonical age period cohort parametrisation is used.  This involves double differences of the
time effects. 
The first row of plots are double differences of the time effects.
The next two rows of plots illustrate the representation theorem depending on the choice of \code{type}.
In both cases the sum of the plots add up to the predictor.	 
	\describe{
		\item{"detrend"}{The last row of plots are double sums of double differences detrend so that that each series starts in
			zero and ends in zero.  The corresponding level and (up to) two linear trends are shown in the middle row of plots.
			The linear trends are identified to be 0 for age, period or cohort equal to its smallest value.
			See note 2 below.
			}
		\item{"sum.sum"}{The last row of plots are double sums of double differences anchored as in the derivation of
			Nielsen (2014b).
			The corresponding level and (up to) two linear trends are shown in the middle row of plots.
			The linear trends are identified to be 0 for the anchoring point U of age, period or cohort as
			described in 
			Nielsen (2014b).
			See note 1 below.
			}
	}

\bold{Model designs with 2 factors.}
If \code{type} is "dif" the canonical two factor parametrisation is used.
This involves single differences.
It is only implemented for \code{model.design} of "AC", "AP", "PC".
It does not apply for \code{model.design} of "APC" because single differences are not identified.
It does not apply for the drift models where \code{model.design} is "Ad", "Pd", "Cd", "t" because it is not clear which time scale the second linear trend should be attributed to.
It is not implemented for \code{model.design} of "tA, "tP", "tC", "1".
The first row of plots are single differences of the time effects.
The next two rows of plots illustrate the representation theorem.  In the second row the level is given and in
the third row plots of single sums of single differences are given, normalised to start in zero.

\bold{Appearance may vary.}
Note, the plots "detrend" and "dif" can give very different appearance of the time effects. The "dif" plots are dominated by
linear trends.  They can therefore be more difficult to interpret than the "detrend" plots, where linear trends are set aside.

\bold{Standard deviations.}
All plots include plots of 1 and 2 standard deviations.  The only exception is the intercept in the case 
\code{model.family} is "poisson.response" as this uses a multinomial sampling scheme, where the intercept is set to increase
in the asymptotic experiment.  The default is to plot standard deviations around zero, so that they represent
a test for zero values of the parameters.
Using the argument \code{sdv.at.zero} the standard deviations can be centered around the estimates. This can give a
very complicated appearance.

\bold{Values of coefficients.}
These can be found using
\code{\link{apc.identify}}.
}
\usage{apc.plot.fit(apc.fit.model,scale=FALSE,
					sdv.at.zero=TRUE,type="detrend",
					sub.plot=NULL,main.outer=NULL,main.sub=NULL,
					cex=NULL,cex.axis=NULL)}
\arguments{
  \item{apc.fit.model}{List. See \code{\link{apc.fit.model}} for a description of the format.}
  \item{scale}{Optional. Logical. If (TRUE) FALSE use scale of (inverse) link function.  Default is FALSE.}
  \item{sdv.at.zero}{Optional. Logical. If FALSE/TRUE standard deviations are plotted around estimates/zero.  Default is TRUE.}
  \item{type}{Optional.  Character.  If "detrend" double sums start and end in zero.  If "sum.sum" double sums anchored as discussed in Nielsen (??). Default is "detrend".}
  \item{sub.plot}{Optional.  Character: "a","b",...,"i".  Only the indicated sub plot is plotted. Default is NULL so all plots shown.}
  \item{main.outer}{Optional.  Character.  Main title in outer margin.  Default is generated internally.}
  \item{main.sub}{Optional.  Vector of 9 characters.  Main titles for individual plots. Default is generated internally.}
  \item{cex}{Optional.  Plot parameter, see \code{\link{par}}. Controls size of text. Default is 1.}
  \item{cex.axis}{Optional.  Plot parameter, see \code{\link{par}}. Controls size of axis annotations. Default is 1.}
}
\note{

(1)
The type "sum.sum" (same as "ss.dd") gives double sums anchored 
to be zero in the three points where
\code{age=cohort=U},
\code{age=U+1,cohort=U}
\code{age=U,cohort=U+1}
with
\code{apc.fit.model$U}
and where
\code{U}
is the integer value of
\code{(per.zero+3)/2}
This corresponds to the representation in
Nielsen (2014b).
The linear plane is parametrised in terms of 
a level, which is the value of the predictor at
\code{age=cohort=U};
an age slope, which is the difference of the values of the predictor at
\code{age=U+1,cohort=U}
and 
\code{age=cohort=U};
an cohort slope, which is the difference of the values of the predictor at
\code{age=U,cohort=U+1}
and 
\code{age=cohort=U}.

(2)
The type "detrend" gives double sums that start in zero and end in zero.
The linear plane is parametrised in terms of
a level, which is the value of the predictor at 
\code{age=cohort=1}, which is usually outside the index set for the data;
while age and cohort slopes are adjusted for the ad hoc identification of the time effects.
}
\references{
Kuang, D., Nielsen, B. and Nielsen, J.P. (2008) Identification of the age-period-cohort model and the extended chain ladder model. Biometrika 95, 979-986. \emph{Download}: \href{http://biomet.oxfordjournals.org/cgi/reprint/95/4/979}{Article}; Earlier version \href{http://www.nuffield.ox.ac.uk/economics/papers/2007/w5/KuangNielsenNielsen07.pdf}{Nuffield DP}.

Nielsen, B. (2014b) Deviance analysis of age-period-cohort models. Work in progress.
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 12 Apr 2015}
\seealso{
\code{\link{data.asbestos}}
and
\code{\link{data.Italian.bladder.cancer}} for information on the data used in the example.

Values of coefficients can be found using \code{\link{apc.identify}}.

The \href{../doc/index.html}{vignette} \href{../doc/Identification.pdf}{Identification.pdf}
has information on the identification.
}
\examples{
#####################
#	Example with Italian bladder cancer data
#	Note that the model.design "AC" cannot be rejected against "APC"
#		so there is little difference between the two plots of those fits.

data.list	<- data.Italian.bladder.cancer()
apc.fit.table(data.list,"poisson.dose.response")
fit.apc		<- apc.fit.model(data.list,"poisson.dose.response","APC")
apc.plot.fit(fit.apc)
dev.new()
fit.ac		<- apc.fit.model(data.list,"poisson.dose.response","AC")
apc.plot.fit(fit.ac)

#	to check the numerical values for the last two rows of plots use
apc.identify(fit.ac)$coefficients.detrend
apc.identify(fit.ac)$coefficients.detrend

#	to get only a sub plot and playing with titles
apc.plot.fit(fit.ac,sub.plot="a",main.outer="My outer title",main.sub="My sub title")
#	to get only a all plots and playing with titles
apc.plot.fit(fit.ac,main.outer="My outer title",main.sub=c("1","2","3","4","5","6","7","8","9"))
}