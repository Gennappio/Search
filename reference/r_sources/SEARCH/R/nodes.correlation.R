nodes.correlation <- function(seeds, x, Blist, covariates = NULL, correlation.threshold = 0.5, correlation.method = "spearman", correlation.type = "standard", correlation.use = "everything") {

	# If doing spearman correlation, we want to turn data into ranks
	#if(correlation.method == "spearman" ) {
	#	x <- t(apply(x, 1, rank));
	#}

	if(correlation.type == "standard") {

		# toConsider <- x[rownames(x) %in% Blist, ];
		# thisProbeset <- x[seed, ];

		# run down the rows of the matrix and get the correlations with x
		#all <- apply(
		#	toConsider,
		#	1,
		#	function(p,q) {
		#		cor(q, p, method = correlation.method, use = "pairwise.complete.obs")
		#	},
		#	q = thisProbeset
		#);

		# filtered <- all[abs(all) > correlation.threshold];

		x.slimmed <- t(x[rownames(x) %in% Blist, ]);
		# use = "pairwise.complete.obs" makes it extremely^100 slow, and hence removed
		x.cor <- cor(x.slimmed[, seeds], x.slimmed, method = correlation.method, use = correlation.use);
		x.cor.list <- as.list(data.frame(t(x.cor)));
		x.cor.list <- lapply(x.cor.list, function(p) { names(p) <- colnames(x.cor); return (p); })

		x.cor.filtered <- lapply(x.cor.list, function(q) { q[abs(q) > correlation.threshold] } )
		names(x.cor.filtered) <- seeds;
		cat(".");

# cat("\nGET CORRELATED - CORRELATION:\n");
# print(str(filtered));
# print (filtered);

# cat("\nIDENTICAL");
# print(identical(gsub("X", "", names(x.cor.filtered)), seeds));
# print(x.cor.filtered);

		# return(filtered);
		return(
			list(
				"CL" = x.cor.filtered, 
				"other.data" = NULL
				)
			);
		}

	stop("correlation.type not supported!");
}
