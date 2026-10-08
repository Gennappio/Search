\name{get.adjacency.matrices}
\alias{get.adjacency.matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.adjacency.matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.adjacency.matrices(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data matrix, each row represents a network}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{list of matrices, each representing a row of input matrix}
\references{
%% ~put references to the literature/web site here ~
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\author{Syed Haider}
\examples{
x <- matrix(data = sample(c(1, 1, NA), 32, replace = TRUE), nrow = 4);
rownames(x) <- c("TP53", "BRCA1", "AURKA", "MYC");
colnames(x) <- c("TP53", "BRCA1", "AURKA", "MYC", "RB", "SNORD4", "ITGA5", "ALK");
adjacency.matrices <- get.adjacency.matrices(x);
}

