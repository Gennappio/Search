get.adjacency.matrices <- function(x) {

	adjacency.matrices.list <- list();
	seeds <- rownames(x);
	seeds.i <- 1;

	# convert every row (subnetwork) to adjacency matrix
	apply(
		x,
		1,
		FUN  = function(y) {
			features <- union(seeds[seeds.i], names(which(!is.na(y))));
			adjacency.matrix <- matrix(
				data = 0, 
				nrow = length(features), 
				ncol = length(features), 
				dimnames = list(features, features)
				);

			adjacency.matrix[seeds[seeds.i], features] <- y[features];
			adjacency.matrix[features, seeds[seeds.i]] <- y[features];

			adjacency.matrices.list[[seeds[seeds.i]]] <<- adjacency.matrix;
			seeds.i <<- seeds.i + 1;
			}
		);

	return (adjacency.matrices.list);
}