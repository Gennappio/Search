write.data.files <- function(results.obj = NULL, data.directory = ".", iteration = 1) {

	# write correlation data of subnetworks with the correlations in a matrix
	for (iter.i in names(results.obj[["correlationMatrix"]])) {

		results.obj[["correlationMatrix"]][[as.character(iter.i)]] <- matrix(
			sprintf("%.3f", results.obj[["correlationMatrix"]][[as.character(iter.i)]]), 
			nrow = nrow(results.obj[["correlationMatrix"]][[as.character(iter.i)]]), 
			dimnames = list(
				rownames(results.obj[["correlationMatrix"]][[as.character(iter.i)]]),
				colnames(results.obj[["correlationMatrix"]][[as.character(iter.i)]])
				)
			);

		# convert "NA" to NA
		tmp.rownames <- rownames(results.obj[["correlationMatrix"]][[as.character(iter.i)]]);
		tmp.colnames <- colnames(results.obj[["correlationMatrix"]][[as.character(iter.i)]]);
		results.obj[["correlationMatrix"]][[as.character(iter.i)]][
			results.obj[["correlationMatrix"]][[as.character(iter.i)]] == "NA"
			] <- NA;

		# transpose and performing row-wise apply to ensure 1xM matrix remains a matrix
		results.obj[["correlationMatrix"]][[as.character(iter.i)]] <- t(
			apply(
				results.obj[["correlationMatrix"]][[as.character(iter.i)]], 
				1, 
				as.numeric
				)
			);
		rownames(results.obj[["correlationMatrix"]][[as.character(iter.i)]]) <- tmp.rownames;
		colnames(results.obj[["correlationMatrix"]][[as.character(iter.i)]]) <- tmp.colnames;

		write.table(
			x = results.obj[["correlationMatrix"]][[as.character(iter.i)]],
			file = paste(data.directory, "subnetworks_correlation__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
			row.names = TRUE,
			col.names = NA,
			quote = FALSE,
			sep = "\t"
			);

		# compress
		system( paste("gzip --force --best ", data.directory, "subnetworks_correlation__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
		}

	# write the shared neighbourhood measures
	for (iter.i in names(results.obj[["O"]])) {

		results.obj[["O"]][[as.character(iter.i)]] <- matrix(
			sprintf("%.3f", results.obj[["O"]][[as.character(iter.i)]]), 
			nrow = nrow(results.obj[["O"]][[as.character(iter.i)]]), 
			dimnames = list(
				rownames(results.obj[["O"]][[as.character(iter.i)]]), 
				colnames(results.obj[["O"]][[as.character(iter.i)]])
				)
			);

		# convert "NA" to NA
		tmp.rownames <- rownames(results.obj[["O"]][[as.character(iter.i)]]);
		tmp.colnames <- colnames(results.obj[["O"]][[as.character(iter.i)]]);

		# transpose and performing row-wise apply to ensure 1xM matrix remains a matrix
		results.obj[["O"]][[as.character(iter.i)]] <- t(
			apply(
				results.obj[["O"]][[as.character(iter.i)]], 
				1, 
				as.numeric
				)
			);
		rownames(results.obj[["O"]][[as.character(iter.i)]]) <- tmp.rownames;
		colnames(results.obj[["O"]][[as.character(iter.i)]]) <- tmp.colnames;

		write.table(
			x = results.obj[["O"]][[as.character(iter.i)]],
			file = paste(data.directory, "shared_neighbourhood__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
			row.names = TRUE,
			col.names = NA,
			quote = FALSE,
			sep = "\t"
			);

		# compress
		system( paste("gzip --force --best ", data.directory, "shared_neighbourhood__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
		}

	# write D (network connectivity) files to file system
	for (iter.i in names(results.obj[["D"]])) {

		results.obj[["D"]][[as.character(iter.i)]] <- sort(
			round(
				results.obj[["D"]][[as.character(iter.i)]],
				digits = 5
				),
			decreasing = TRUE
			);
		write.table(
			x = results.obj[["D"]][[as.character(iter.i)]],
			file = paste(data.directory, "connectivity__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
			row.names = TRUE,
			col.names = NA,
			quote = FALSE,
			sep = "\t"
			);

		# compress
		system( paste("gzip --force --best ", data.directory, "connectivity__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
		}

	# write the network growth files to the file system
	for (iter.i in names(results.obj[["network.growth"]])) {

		write.table(
			x = results.obj[["network.growth"]][[as.character(iter.i)]],
			file = paste(data.directory, "network_growth__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
			row.names = TRUE,
			col.names = NA,
			quote = FALSE,
			sep = "\t"
			);

		# compress
		system( paste("gzip --force --best ", data.directory, "network_growth__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
		}

	# write the node degree file
	for (iter.i in names(results.obj[["nodes.degree"]])) {

		results.obj[["nodes.degree"]][[as.character(iter.i)]] <- round(
			results.obj[["nodes.degree"]][[as.character(iter.i)]], digits = 5
			);

		write.table(
			x = sort(results.obj[["nodes.degree"]][[as.character(iter.i)]], decreasing = TRUE),
			file = paste(data.directory, "nodes_degree__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
			row.names = TRUE,
			col.names = NA,
			quote = FALSE,
			sep = "\t"
			);

		# compress
		system( paste("gzip --force --best ", data.directory, "nodes_degree__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
		}

	# write the node degree file
	if (!is.null(results.obj[["other.data"]])) {
		for (iter.i in names(results.obj[["other.data"]])) {
			write.table(
				x = results.obj[["other.data"]][[as.character(iter.i)]],
				file = paste(data.directory, "other_data__iteration_", iteration, "_", iter.i, ".txt", sep = ""),
				row.names = TRUE,
				col.names = NA,
				quote = FALSE,
				sep = "\t"
				);

			# compress
			system( paste("gzip --force --best ", data.directory, "other_data__iteration_", iteration, "_", iter.i, ".txt", sep = "") );
			}
		}

	return (results.obj);
}
