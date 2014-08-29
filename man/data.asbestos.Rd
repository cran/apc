\name{data.asbestos}
\alias{data.asbestos}
\title{Asbestos data}
\description{
Function that organises asbestos data in \code{\link{apc.data.list}} format.

Counts of mesothelioma deaths in the UK by age and period.   Mesothelioma is most often caused by exposure to asbestos.

The data set is in "PA"-format.

The primary data set includes ages 25-89, which is obtained when using the function without arguments or with argument \code{all.age.groups=FALSE}. The secondary data includes younger and older age groups, which is obtained when using the function with argument \code{all.age.groups=TRUE}. The \code{apc.package} is at present not aimed at such unbalanced data.
}
\usage{data.asbestos(all.age.groups = FALSE)}
\arguments{
  	\item{all.age.groups}{logical.  If FALSE (default), only age groups 25-89 are included.}
}
\value{The value is a list in \code{\link{apc.data.list}} format.
	\item{response}{matrix of cases. Numbers of mesothelioma deaths by period and age. Period runs 1967-2007. Age runs 25-89 when \code{all.age.groups=FALSE}.  "PA"-format.}
	\item{dose}{NULL}
	\item{data.format}{logical equal to "PA".  Data organised with period-groups in rows and age-groups in columns.}
	\item{age1}{numeric equal to 25.  This is the label for the first age group of 25.}
	\item{per1}{numeric equal to 1967.  This is the label for the first period group of 1967.}
	\item{coh1}{NULL.  Not needed when data.format="PA"}
	\item{unit}{numeric equal to 1.  This is the width of the age and period groups.}
	\item{per.zero}{NULL.  Not needed when data.format="PA"} 
	\item{per.max}{NULL.  Not needed when data.format="PA"}
	\item{time.adjust}{0.  Thus age=89 in period=1967 corresponds to cohort=1967-89+0=1878.}
}
\source{Data were prepared for the Asbestos Working Party by the UK Health and Safety Executive.
An APC analysis of these data can be found in Martinez Miranda, Nielsen and Nielsen (2014)}
\references{Martinez Miranda, M.D., Nielsen, B. and Nielsen, J.P. (2014)
Inference and forecasting in the age-period-cohort model with unknown
exposure with an application to mesothelioma mortality.
To appear in \emph{Journal of the Royal Statistical Society series} A.}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 29 August 2014}
\seealso{General description of \code{\link{apc.data.list}} format.
}
\examples{
#########################
#	apc data list

data.list	<- data.asbestos()
objects(data.list)

#####################
#	Figure 1,a-c from
#	Miranda Martinex, Nielsen and Nielsen (2013).

data.list	<- data.asbestos()	
apc.plot.data.sums(data.list,type="l")

#####################
#	Figure 1,d from
#	Miranda Martinex, Nielsen and Nielsen (2013).
data.list	<- data.asbestos()	
apc.plot.data.within(data.list,type="l",lty=1)

}