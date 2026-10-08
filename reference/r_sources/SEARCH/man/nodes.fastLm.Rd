\name{nodes.fastLm}
\alias{nodes.fastLm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nodes.fastLm
}
\description{
C implementation of linear models - much faster compared to glm()
}
\usage{
nodes.fastLm(seeds, x, Blist, covariates = NULL, significance.threshold = 0.01, adjust.p = NULL, glm.family = "gaussian", formula.string = NULL, mc.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seeds}{seeds set}
  \item{x}{data matrix}
  \item{Blist}{universe set}
  \item{covariates}{covariates matrix to be used for additional model components in formula.string. The order of rows should be strictly the same as order of columns of data in SEARCH.randomise(...)}
  \item{significance.threshold}{P value threshold to be applied for choosing that seeds significantly associated with other nodes}
  \item{adjust.p}{specify method name to use to correct P values in each iteration for multiple comparisons. defaults to no correction. values are same as accepted ?p.adjust}
  \item{glm.family}{distribution to be used by glm. Only gaussian is supported}
  \item{formula.string}{contains the additional model variables in addition to the base model y ~ x.node (as in R)\cr eg. a model y ~ x.node + v + x.node:v should be specified as "v + x.node:v"}
  \item{mc.cores}{number of multiple cores to be used for parallelisation. NOTE: Windows does not support this function. OS X is often faced by a warning Tcl_WaitForEvent: Notifier not initialized, in which case, stick to mc.core = 1 which is default}
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
