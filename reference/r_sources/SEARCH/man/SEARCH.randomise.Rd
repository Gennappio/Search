\name{SEARCH.randomise}
\alias{SEARCH.randomise}
\title{
SEARCH.randomise
}
\description{
Wrapper method to call seed clustering algorithm
}
\usage{
SEARCH.randomise(data, seed.nodes, universe.nodes, covariates = NULL, scale = TRUE, D.lower = 0.1, D.higher = 0.24, O.lower = 0.01, keep.seeds = TRUE, max.iter = 1, verbose = TRUE, sd.threshold = 0.1, n.boot = 10, n.boot.start = 1, sampling = "bootstrap", sampling.seed = NULL, output.directory = ".", nodes.association.function = "nodes.correlation", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix}
  \item{seed.nodes}{seeds set}
  \item{universe.nodes}{universe set}
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
  \item{keep.seeds}{If FALSE, apply O.lower based filtering on all seeds, new and old. Otherwise, only on new seeds}
  \item{max.iter}{Integer indicating how many segments to travel away from seeds for new recruits}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{sd.threshold}{
%%     ~~Describe \code{sd.threshold} here~~
}
  \item{n.boot}{
%%     ~~Describe \code{n.boot} here~~
}
  \item{n.boot.start}{An integer to specify where to start the iteration counter. Defaults to 1}
  \item{sampling}{A string specifying name of sampling method}
  \item{sampling.seed}{deterministic/reproducible sampling. Critical for matched sampling across tumour types/subtypes}
  \item{output.directory}{Path of the output directory where results will be stored}
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

# get data directory 
data.directory <- get.program.defaults()[["test.data.dir"]];
output.directory.corr <- "./seedclusteringoutputCORR/";
output.directory.glm <- "./seedclusteringoutputGLM/";

# initialise params
feature.selection.dataset <- c("Breastdata1");

# read in expression data
x <- as.matrix(
  read.table(
    file = paste(data.directory, feature.selection.dataset, "/mRNA_abundance.txt.gz" , sep = ""),
    row.names = 1,
    header = TRUE,
    sep = "\t"
    )
  );

# choose seed nodes
seed.nodes <- rownames(x)[1:10];

# choose universe nodes
universe.nodes <- rownames(x);

# run interface method (correlation based)
results.obj <- SEARCH.randomise(
  data = x,
  seed.nodes = seed.nodes,
  universe.nodes = universe.nodes,
  covariates = NULL,
  scale = TRUE, 
  D.lower = 0, 
  D.higher = 0, 
  O.lower = 0.1,
  keep.seeds = TRUE,
  max.iter = 3,
  verbose = TRUE, 
  sd.threshold = 0,
  n.boot = 2,
  sampling.seed = 1,
  output.directory = output.directory.corr,
  nodes.association.function = "nodes.correlation",
  correlation.threshold = 0.5,
  correlation.method = "spearman",
  correlation.type = "standard"
  );

# run interface method (GLM based)
results.obj <- SEARCH.randomise(
  data = x,
  seed.nodes = seed.nodes,
  universe.nodes = universe.nodes,
  covariates = NULL,
  scale = TRUE, 
  D.lower = 0, 
  D.higher = 0, 
  O.lower = 0.1,
  keep.seeds = TRUE,
  max.iter = 2,
  verbose = TRUE, 
  sd.threshold = 0,
  n.boot = 2,
  sampling.seed = 1,
  output.directory = output.directory.glm,
  nodes.association.function = "nodes.glm",
  significance.threshold = 1e-05,
  glm.family = "gaussian",
  mc.cores = 5
  );

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
\keyword{ networks }
