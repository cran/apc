\name{apc.data.list}
\alias{apc.data.list}
\title{Arrange data as an apc.data.list}
\description{This is step 1 of the apc analysis.

The apc package is aimed at range of data types.  This analysis and labelling of parameters depends on
the choice data type.  In order to keep track of this choice the data first has to be arranged as an
apc.data.list.  The function purpose of this function is to aid the user in constructing a list with the right information.

Age period cohort analysis is used in two situations.  A dose-response situation, where both doses (exposure, risk set, cases) and
responses (counts of deaths, outcomes) are available.  And a response situation where only a response is available.  If
the aim is to directly model mortality ratios (counts of death divided by exposure) this will be thought of a response

The \code{apc.data.list} gives sufficient information for the further analysis. It is sufficient to store this information.
It has 2 obligatory arguments, which are a response matrix and a character indicating the data format.
It also has some further optional arguments, which have certain default values. 
Some times it may be convenient to add further arguments to the \code{apc.data.list}.  This will not affect the apc analysis.
}
\usage{apc.data.list(response, data.format, dose=NULL,
					age1=1, per1=1, coh1=1, unit=1,
					per.zero=NULL, per.max=NULL, time.adjust=0)}
\arguments{
	\item{response}{matrix (or vector).  Numbers of responses.  It should have a format matching \code{data.format}}
	\item{data.format}{character.  The following options are implemented:
  	\describe{
			\item{"AC"}{has    age/cohort as increasing row/column index.}
			\item{"AP"}{has    age/period as increasing row/column index.}
			\item{"CA"}{has cohort/age    as increasing row/column index.}
			\item{"CL"}{has cohort/age 	  as increasing row/column index, triangular.}
			\item{"CP"}{has cohort/period as increasing row/column index.}
			\item{"PA"}{has period/age    as increasing row/column index.}
			\item{"PC"}{has period/cohort as increasing row/column index.}
			\item{"trapezoid"}{has    age/period as increasing row/column index, period-diagonals are NA for period <= per.zero and >per.zero+per.max.}
		}}
 	\item{dose}{\emph{Optional}. matrix or NULL.  Numbers of doses.  It should have same format as \code{response}.}
    \item{age1}{\emph{Optional}. Numeric or NULL.  Time label for youngest age group.  Used if \code{data.format} is "AC", "AP", "CA", "CL", "CL.vector.by.row", "PA", "trapezoid". If NULL default is 1.}
    \item{per1}{\emph{Optional}. Numeric or NULL.  Time label for oldest period group.  Used if \code{data.format} is "AP", "CP", "PA", "PC". If NULL default is 1.}
    \item{coh1}{\emph{Optional}. Numeric or NULL.  Time label for youngest age group.  Used if \code{data.format} is "AC", "CA", "CL", "CL.vector.by.row", "CP", "PC", "trapezoid". If NULL default is 1.}
    \item{unit}{\emph{Optional}. Numeric or NULL.  Common time steps for age, period and cohort.  If NULL default is 1.}
	\item{per.zero}{\emph{Optional}. Numeric or NULL.  Needed if data format is "trapezoid".}
	\item{per.max}{\emph{Optional}. Numeric or NULL.  Needed if data format is "trapezoid".}
	\item{time.adjust}{\emph{Optional}. Numeric.  Time labels are based on two of age1, per1 and coh1.  The third time label is computed according to the formula age1+coh1=per1+time.adjust.  Default is 0.  If age1=coh=1 it is natural to choose time.adjust=1.}
}
\value{
	  	\item{response}{matrix (or vector).  Numbers of responses.}
	 	\item{dose}{matrix (or NULL).  Numbers of doses.}
		\item{data.format}{character.}
	    \item{age1}{Numeric. Default of 1.}		
	    \item{per1}{Numeric. Default of 1.}		
	    \item{coh1}{Numeric. Default of 1.}		
	    \item{unit}{Numeric.  Default of 1.}		
		\item{per.zero}{Numeric.  If data.format is not "trapezoid" the value is NULL. If data.format is "trapezoid" the default is \code{per.zero=0}.}
		\item{per.max}{Numeric.  If data.format is not "trapezoid" the value is NULL. If data.format is "trapezoid" the default is \code{per.max=nrow(response)+ncol(response)-1-per.zero}.}
		\item{time.adjust}{Numeric.  Default of 0.}
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 18 Mar 2015}
\seealso{The below example shows how the \code{\link{data.Japanese.breast.cancer}} data.list was generated.
Other provided data sets include
\code{\link{data.asbestos}}
\code{\link{data.Belgian.lung.cancer}}
\code{\link{data.Italian.bladder.cancer}}.

A subset of the data can be selected using \code{\link{apc.data.list.subset}}.
}
\examples{
###############
#	Artificial data
#	(1) Generate a 5x7 matrix and make arbitrary decisions for rest

response <- matrix(data=seq(1:35),nrow=5,ncol=7)
data.list	<- list(response=response,dose=NULL,data.format="AP",
					age1=25,per1=1955,coh1=NULL,unit=5,
					per.zero=NULL,per.max=NULL)
data.list

#	(2) Chain Ladder data

k			<- 5
v.response 	<- seq(1:(k*(k+1)/2))
data.list	<- apc.data.list(response=vector.2.triangle(v.response,k),
							data.format="CL.vector.by.row",age1=2001)
data.list


###############
#	Japanese breast cancer
#	This is the code used to generate the data.Japanese.breast.cancer
v.rates		<- c( 0.44, 0.38, 0.46, 0.55, 0.68,
			 	  1.69, 1.69, 1.75, 2.31, 2.52,
				  4.01, 3.90, 4.11, 4.44, 4.80,
				  6.59, 6.57, 6.81, 7.79, 8.27,
				  8.51, 9.61, 9.96,11.68,12.51,
				 10.49,10.80,12.36,14.59,16.56,
				 11.36,11.51,12.98,14.97,17.79,
				 12.03,10.67,12.67,14.46,16.42,
				 12.55,12.03,12.10,13.81,16.46,
				 15.81,13.87,12.65,14.00,15.60,
				 17.97,15.62,15.83,15.71,16.52)
v.cases		<- c(   88,   78,  101,  127,  179,
				   299,  330,  363,  509,  588,
				   596,  680,  798,  923, 1056,
				   874,  962, 1171, 1497, 1716,
				  1022, 1247, 1429, 1987, 2398,
				  1035, 1258, 1560, 2079, 2794,
				   970, 1087, 1446, 1828, 2465,
				   820,  861, 1126, 1549, 1962,
				   678,  738,  878, 1140, 1683,
				   640,  628,  656,  900, 1162,
				   497,  463,  536,  644,  865)				 
row.names	<- c("25-29","30-34","35-39","40-44","45-49",
				 "50-54","55-59","60-64","65-69","70-74","75-79")
col.names	<- c("1955-1959","1960-1964","1965-1969","1970-1974","1975-1979")

rates	<- matrix(data=v.rates,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
cases	<- matrix(data=v.cases,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))

# 	A data list is now constructed as follows
#	note that list entry rates is redundant,
#	but included since it represents original data

data.Japanese.breast.cancer	<- list(rates=rates,response=cases,
			dose=cases/rates,data.format="AP",
			age1=25,per1=1955,coh1=NULL,unit=5,
			per.zero=NULL,per.max=NULL,time.adjust=0)
}