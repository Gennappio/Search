get.matrix.format <- function(x, all.features) {

	a <- names(x);
	O <- matrix(data = 0, nrow = length(a), ncol = length(a));
	rownames(O) <- a;
	colnames(O) <- a;

	correlation.matrix <- matrix(
		data = NA, 
		nrow = length(a), 
		ncol = length(all.features),
		dimnames = list(a, all.features)
		);

	for(i in 1:length(a)) {
		correlation.matrix[ a[i], names(x[[a[i]]]) ] <- x[[a[i]]];
		for(j in i:length(a)) {
			O[i,j] <- getO(x[[i]], x[[j]]);
			}
		}

	# remove columns from subnetworks corrrelation matrix which represents features that 
	# are not involved in any subnetwork
	# convert to data.frame just to handle 1xM matrix otherwise apply vectorise it
	correlation.matrix <- as.data.frame(correlation.matrix);
	correlation.matrix <- as.matrix(
		correlation.matrix[, apply(correlation.matrix, 2, function(x) !all(is.na(x)))]
		);

	return( list("O" = O, "subnetworks.correlation" =  correlation.matrix) );
	}
