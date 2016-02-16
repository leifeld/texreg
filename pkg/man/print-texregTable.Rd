\name{print.texregTable}
\alias{print.texregTable}
\alias{texregTable}
\title{Print the output of a screenreg, texreg, or htmlreg call}
\description{Print the output of a screenreg, texreg, or htmlreg call.}
\usage{
\S3method{print}{texregTable}(x, ...)
}
\details{
This function prints a texregTable object, which results from a screenreg, texreg, or htmlreg call. Most of the time, this function is called implicitly by just entering the name of the object.
}
\arguments{
\item{x}{ A texregTable object. This is basically a simple \code{character} object with an additional class name called \code{texregTable}. }
\item{...}{ Additional arguments to be handed over to the \code{cat} function. }
}
\seealso{
\code{\link{texreg-package} \link{texreg}}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.com})
}
