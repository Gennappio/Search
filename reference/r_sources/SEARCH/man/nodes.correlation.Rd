\name{nodes.correlation}
\alias{nodes.correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nodes.correlation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nodes.correlation(seeds, x, Blist, covariates = NULL, correlation.threshold = 0.5, correlation.method = "spearman", correlation.type = "standard", correlation.use = "everything")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seeds}{seeds set}
  \item{x}{data matrix}
  \item{Blist}{universe set}
  \item{covariates}{covariates matrix to be used for additional model components in formula.string. The order of rows should be strictly the same as order of columns of data in SEARCH.randomise(...)}
  \item{correlation.threshold}{
%%     ~~Describe \code{t} here~~
}
  \item{correlation.method}{
%%     ~~Describe \code{correlation.method} here~~
}
  \item{correlation.type}{
%%     ~~Describe \code{correlation.type} here~~
}
  \item{correlation.use}{use param of cor() function in R. Note that "pairwise.complete.obs" is extremely slow}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Syed Haider, Francesca M. Buffa
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
