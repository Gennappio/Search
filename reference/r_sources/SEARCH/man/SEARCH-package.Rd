\name{SEARCH-package}
\alias{SEARCH-package}
\docType{package}
\title{
Seed Agglomerative and Recursive Clustering with Hypothesis Oriented Initialization (SEARCH)
}
\description{
Creates clusters of subnetwork features around pre-specified seeds

}
\details{
\tabular{ll}{
Package: \tab SEARCH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-10\cr
License: \tab GPL-2\cr
}

}

\author{
Syed Haider, Francesca M. Buffa
Maintainer: Syed Haider, Francesca M. Buffa\cr
MRC Weatherall Institute of Molecular Medicine\cr
Department of Oncology\cr
University of Oxford.
}
\references{
Buffa FM et al. (2010) Large meta-analysis of multiple cancers reveals a common, compact and highly prognostic hypoxia metagene. British Journal of Cancer
}
\keyword{ clustering }
\seealso{
}
\examples{

# get data directory 
data.directory <- get.program.defaults()[["test.data.dir"]];
output.directory.corr <- "./seedclusteringoutputCORR/";
output.directory.glm <- "./seedclusteringoutputGLM/";
output.directory.glm.covariates <- "./seedclusteringoutputGLMCovariates/";
output.directory.fastLm.covariates <- "./seedclusteringoutputFastLMCovariates/";

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

# run interface method (GLM based with covariates)
covariates <- data.frame("ER" = (1:ncol(x)) \%\% 2);
results.obj <- SEARCH.randomise(
  data = x,
  seed.nodes = seed.nodes,
  universe.nodes = universe.nodes,
  covariates = covariates,
  formula.string = "ER",
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
  output.directory = output.directory.glm.covariates,
  nodes.association.function = "nodes.glm",
  significance.threshold = 1e-05,
  glm.family = "gaussian",
  mc.cores = 5
  );

# run interface method (Fast linear models with covariates)
covariates <- data.frame("ER" = (1:ncol(x)) \%\% 2);
results.obj <- SEARCH.randomise(
  data = x,
  seed.nodes = seed.nodes,
  universe.nodes = universe.nodes,
  covariates = covariates,
  formula.string = "ER",
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
  output.directory = output.directory.fastLm.covariates,
  nodes.association.function = "nodes.fastLm",
  significance.threshold = 1e-05,
  glm.family = "gaussian",
  mc.cores = 5
  );

# (optional) create network plots

}
