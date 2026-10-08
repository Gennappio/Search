\name{make.networkplot}
\alias{make.networkplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.networkplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.networkplot(adjacency.matrix, network.layout = "layout.fruchterman.reingold", file.name = NULL, height = 5, width = 5, main.title = NULL, resolution = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacency.matrix}{
%%     ~~Describe \code{adjacency.matrix} here~~
}
  \item{network.layout}{A string representing layout as supported in igraph layout, defaults to layout.fruchterman.reingold}
  \item{file.name}{
%%     ~~Describe \code{file.name} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{main.title}{
%%     ~~Describe \code{main.title} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
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
\author{Syed Haider}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
output.directory <- "./seedclusteringoutput/graphs/";

# prepare random dataset
x <- matrix(data = sample(c(1, -1, 0), 32, replace = TRUE), nrow = 4);
rownames(x) <- c("TP53", "BRCA1", "AURKA", "MYC");
colnames(x) <- c("TP53", "BRCA1", "AURKA", "MYC", "RB", "SNORD4", "ITGA5", "ALK");
adjacency.matrices <- get.adjacency.matrices(x);

# create network 
make.networkplot(
  adjacency.matrix = adjacency.matrices[[1]], 
  file.name = paste(output.directory, "cancer_genes_1.tiff", sep = "") 
  );

# create network (circle layout)
make.networkplot(
  adjacency.matrix = adjacency.matrices[[2]], 
  network.layout = "layout.circle",
  file.name = paste(output.directory, "cancer_genes_2.tiff", sep = "") 
  );
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ networks }
\keyword{ graphs }
