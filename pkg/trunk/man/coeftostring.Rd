\name{coeftostring}
\alias{coeftostring}
\title{Convert coefficients into standardized strings}
\description{Convert coefficients into standardized strings.}
\usage{
coeftostring(x, lead.zero = FALSE, digits = 2)
}
\details{
This function converts numbers into strings and standardizes 
them according to some simple rules. The function is used by 
the texreg function.
}
\arguments{
\item{x}{ A numeric object, for example a coefficient resulting from a regression model. }
\item{lead.zero}{ If the number starts with "0" or "-0": should the zero be retained or removed? If true, the leading zero is kept. }
\item{digits}{ The number of decimal places to be used. }
}
\seealso{
\code{\link{texreg-package} \link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})
}
\examples{
coeftostring(-0.345, lead.zero = FALSE)  # this should give "-.34"
}
