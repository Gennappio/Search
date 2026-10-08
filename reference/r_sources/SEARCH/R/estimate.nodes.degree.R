estimate.nodes.degree <- function(correlation.matrix.list = NULL) {

	nodes.degree.list <- list();

	for (iter.i in names(correlation.matrix.list)) {

		correlation.matrix <- correlation.matrix.list[[as.character(iter.i)]];

		# sanity check
		if (is.null(correlation.matrix)) {
			stop("\ncant proceed as correlation matrix is NULL")
			}	

		# to avoid where a seed is not selected by another seed. in such cases, it will not be
		# present in the columns of correlationMatrix
		common.nodes <- intersect(rownames(correlation.matrix), colnames(correlation.matrix));
		correlation.matrix <- as.matrix(correlation.matrix[common.nodes, common.nodes]);

		# number of other nodes with strong correlation
		nodes.degree <- apply(correlation.matrix, 1, FUN = function(x) { length(which(x != "NA")) } );

		# fraction of total nodes
		nodes.degree <- nodes.degree/ncol(correlation.matrix);

		# assign names
		names(nodes.degree) <- rownames(correlation.matrix);

		# populate return data structure
		nodes.degree.list[[as.character(iter.i)]] <- nodes.degree;

		}

	return (nodes.degree.list);
	}