\name{get.program.defaults}
\alias{get.program.defaults}
\title{A utility function to return the inst/ directory of the installed package and other default settings}
\description{A utility function to return the inst/ directory of the installed package to get the test datasets and other program related data contents}
\usage{
get.program.defaults(entry.secret = "hypoxia_angiogenesis_sub_networks_all.txt");
}
\arguments{
  \item{entry.secret}{Name of a file that we know exists in the package}
}
\value{Returns path to the directory where the contents of this package are installed and names of other input files}
\author{Syed Haider}
\seealso{}
\examples{
program.data <- get.program.defaults();
}
\keyword{IO}
