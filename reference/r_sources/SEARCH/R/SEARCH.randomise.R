SEARCH.randomise <- function(data, seed.nodes, universe.nodes, covariates = NULL, scale = TRUE, D.lower = 0.1, D.higher = 0.24, O.lower = 0.01, keep.seeds = TRUE, max.iter = 1, verbose = TRUE, sd.threshold = 0.1, n.boot = 10, n.boot.start = 1, sampling = "bootstrap", sampling.seed = NULL, output.directory = ".", nodes.association.function = "nodes.correlation", ...) {

	x <- data;

	# sanity checks
	if(!is.matrix(x) || nrow(x) < 2 || ncol(x) < 2) {
		stop("cant proceed, data is either not a matrix OR does not contain data");
		}

	# restrict seed.nodes to features present in expression data
	seed.nodes <- intersect(seed.nodes, rownames(x));
	if(length(seed.nodes) == 0) {
		stop("cant proceed, seed.nodes and rownames in data have no overlap");
		}

	# warn that data should be on comparable scale when using nodes.glm
	if (nodes.association.function == "nodes.glm") {
		cat("\nWARNING: Please scale data (if you have not done so) when using nodes.glm(..) so that coefficients/betas are comparable to nodes.correlation(..) method if you want to compare these later");
		}

	# run seed clustering algorithm n.boot times
	i <- n.boot.start:(n.boot.start+n.boot-1);
	return.obj <- lapply(
		i, 
		.SEARCH.randomise.inner, 
		x = x, 
		Alist = seed.nodes, 
		Blist = universe.nodes, 
		covariates = covariates,
		scale = scale, 
		D.lower = D.lower, 
		D.higher = D.higher, 
		O.lower = O.lower, 
		keep.Alist = keep.seeds,
		max.iter = max.iter, 
		verbose = verbose,
		sd.threshold = sd.threshold, 
		sampling = sampling,
		sampling.seed =  sampling.seed,
		output.directory = output.directory,
		nodes.association.function = nodes.association.function,
		...
		);

	return(return.obj);
}
