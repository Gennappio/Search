\name{SEARCH}
\alias{SEARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SEARCH
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SEARCH(x, Alist, Blist, covariates = NULL, scale = TRUE, D.lower = 0, D.higher = 0, O.lower = 0, keep.Alist = TRUE, max.iter = 1, sd.threshold = 0, verbose = TRUE, nodes.association.function = "nodes.correlation", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data matrix}
  \item{Alist}{seeds set}
  \item{Blist}{universe set}
  \item{covariates}{covariates matrix to be used for additional model components in formula.string. The order of rows should be strictly the same as order of columns of data in SEARCH.randomise(...)}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{D.lower}{
%%     ~~Describe \code{D.lower} here~~
}
  \item{D.higher}{
%%     ~~Describe \code{D.higher} here~~
}
  \item{O.lower}{
%%     ~~Describe \code{O.lower} here~~
}
  \item{keep.Alist}{If FALSE, apply O.lower based filtering on all seeds, new and old. Otherwise, only on new seeds}
  \item{max.iter}{Integer indicating how many segments to travel away from seeds for new recruits}
  \item{sd.threshold}{
%%     ~~Describe \code{sd.threshold} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{nodes.association.function}{name of the function to be used for estimating network parameters. By default \code{nodes.correlation} is used which estimates correlation between the seed nodes and the universe nodes}
  \item{...}{other params to be passed on to user-defined method for estimating network parameters}
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
\keyword{ clustering }
\keyword{ networks }
