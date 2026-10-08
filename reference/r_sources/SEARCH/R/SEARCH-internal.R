.SEARCH.randomise.inner <- function(iteration = 1, x, Alist, Blist, covariates = NULL, scale, D.lower, D.higher, O.lower, keep.Alist, max.iter, verbose, sd.threshold, sampling, sampling.seed, output.directory = ".", nodes.association.function = "nodes.correlation", ...) {

	# create the output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");
	logs.dir <- paste(output.directory, "/logs/", sep = "");

	if (!file.exists(out.dir)) {
		dir.create(out.dir, recursive = TRUE);
		dir.create(graphs.dir, recursive = TRUE);
		dir.create(logs.dir, recursive = TRUE);
		}

	# set sampling seed if required
	if (!is.null(sampling.seed)) {
		cat("\nsetting sampling seed =", iteration + as.numeric(sampling.seed) - 1);
		set.seed(iteration + as.numeric(sampling.seed) - 1);
		}

	# set sampling type specific parameters
	if (sampling == "bootstrap") {
		idx <- 1:(dim(x)[2]);
		random.samples <- sample(idx, replace = TRUE);
		x <- x[, random.samples];
		covariates.sample <- covariates[random.samples, , drop = FALSE];
		}
	else if(sampling == "jacknife") {
		idx <- 1:(dim(x)[2]);
		random.samples <- sample(idx, size = dim(x)[2]-1, replace = FALSE)
		x <- x[, random.samples];
		covariates.sample <- covariates[random.samples, , drop = FALSE];
		}
	else {
		stop("\nQuitting. unknown sampling method");
		}

	# recover memory
	for (i in 1:5) { gc(); }

	# run seed clustering algorithm
	results.obj <- SEARCH(
		x,
		Alist = Alist,
		Blist = Blist,
		covariates = covariates.sample,
		scale = scale, 
		D.lower = D.lower, 
		D.higher = D.higher, 
		O.lower = O.lower, 
		keep.Alist = keep.Alist,
		max.iter = max.iter, 
		sd.threshold = sd.threshold,
		verbose = verbose,
		nodes.association.function = nodes.association.function,
		...
		);

	# calculate node degree (as a fraction of total nodes in the dataset)
	results.obj[["nodes.degree"]] <- estimate.nodes.degree(
		correlation.matrix.list = results.obj[["correlationMatrix"]]
		);

	# write raw data and default plots to the file system
	results.obj <- write.data.files(results.obj = results.obj, data.directory = out.dir, iteration = iteration);	

	# do minimal plotting
	make.default.plots(results.obj = results.obj, seed.nodes = Alist, graphs.directory = graphs.dir, iteration = iteration);

	# summarise results object and return summary only as everything is written to file system
	summary.obj <- list();
	for (list.element in names(results.obj)) {

		if (list.element %in% c("removed", "correlationLists")) {
			summary.obj[[list.element]] <- length(results.obj[[list.element]]);
			}
		else if (list.element %in% c("correlationMatrix", "O")) {
			summary.obj[[list.element]] <- lapply(results.obj[[list.element]], dim);
			}
		else if (list.element %in% c("D", "network.growth")) {
			summary.obj[[list.element]] <- lapply(results.obj[[list.element]], length);
			}
		else {
			}
		}

	return (summary.obj);
}
