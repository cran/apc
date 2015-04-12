\name{vector.2.triangle}
\alias{vector.2.triangle}
\title{Organise vector as matrix with triangle structure}
\description{
Organise vector as matrix with triangle structure. Useful for reserving data.
It may be easy to input data as vector. This function organises data in "CL"-format,
which is a cohort-age (policy year-development) year matrix
}
\usage{vector.2.triangle(v,k)}
\arguments{
	\item{v}{vector. Length k*(k+1)/2}
	\item{k}{integer. Dimension	    }
}	
\value{The value is a list in \code{\link{apc.data.list}} format.
	\item{m}{matrix with "CL" format Dimension kxk.      
			Upper left triangle filled by v, row by row.
			Remaining entries NA}                        
}
\author{Bent Nielsen <bent.nielsen@nuffield.ox.ac.uk> 7 Feb 2015}
\examples{
#########################

vector.2.triangle(1:10,4)

}
