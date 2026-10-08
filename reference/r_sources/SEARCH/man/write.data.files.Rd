\name{write.data.files}
\alias{write.data.files}
\title{write data object to file system as text files}
\description{write data object to file system as text files}
\usage{
write.data.files(results.obj = NULL, data.directory = ".", iteration = 1);
}
\arguments{
  \item{results.obj}{SEARCH returned results object in list format}
  \item{data.directory}{Path of the output directory where results will be stored}
  \item{iteration}{iteration number of randomisation}
}
\value{\code{results.obj} however in sync with data files written to the file system following rounding to 3-5 decimals places}
\author{Syed Haider}
\examples{
}
\keyword{IO}
