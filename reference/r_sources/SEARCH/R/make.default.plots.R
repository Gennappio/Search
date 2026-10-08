make.default.plots <- function(results.obj = NULL, seed.nodes = NULL, graphs.directory = ".", iteration = 1) {

	# inter-seed correlation network
	# seeds and other nodes correlated with seeds
	for (iter.i in names(results.obj[["correlationMatrix"]])) {

		inter.seed.network <- results.obj[["correlationMatrix"]][[as.character(iter.i)]];

		# to avoid where a seed is not selected by another seed. in such cases, it will not be
		# present in the columns of correlationMatrix
		common.nodes <- intersect(rownames(inter.seed.network), colnames(inter.seed.network));
		inter.seed.network <- as.matrix(inter.seed.network[common.nodes, common.nodes]);

		# ignore self interactions and set NA to 0 bcoz igraphs treat NA as edges in weighted mode
		diag(inter.seed.network) <- 0;
		inter.seed.network[is.na(inter.seed.network)] <- 0;

		# create network 
		if (sum(abs(inter.seed.network)) > 0) {
			make.networkplot(
				adjacency.matrix = inter.seed.network, 
				# network.layout = "layout.circle",
				file.name = paste(graphs.directory, "inter_seed_correlation__iteration_", iteration, "_", iter.i, ".png", sep = ""),
				main.title = paste("Inter seed correlation network (", iter.i,  ")", sep = "")
				);
			}
		else {
			cat("\nSkipping correlation plot as there exists none (all zeros)");
			}
		}

	# recover memory
	for (i in 1:5) { gc(); }

	# inter-seed shared neighbourhood network
	# limit to seed v seed networks only
	for (iter.i in names(results.obj[["O"]])) {

		inter.seed.network <- results.obj[["O"]][[as.character(iter.i)]];
		common.nodes <- intersect(
			seed.nodes,
			intersect(rownames(inter.seed.network), colnames(inter.seed.network))
			);

		inter.seed.network <- as.matrix(inter.seed.network[common.nodes, common.nodes]);

		# ignore self interactions and set NA to 0 bcoz igraphs treat NA as edges in weighted mode
		diag(inter.seed.network) <- 0;
		inter.seed.network[is.na(inter.seed.network)] <- 0;

		# create network 
		if (sum(abs(inter.seed.network)) > 0) {
			make.networkplot(
				adjacency.matrix = inter.seed.network, 
				network.layout = "layout.circle",
				file.name = paste(graphs.directory, "inter_seed_shared_neighbourhood__iteration_", iteration, "_", iter.i, ".png", sep = ""),
				main.title = paste("Inter seed shared neighbourhood network (", iter.i,  ")", sep = "")
				);
			}
		else {
			cat("\nSkipping shared neighbourhood plot as there exists none (all zeros)");
			}
		}

	# create betweeness plot
	for (iter.i in names(results.obj[["nodes.degree"]])) {

		line.plot.data <- list();
		nodes.frequency.table <- table(results.obj[["nodes.degree"]][[as.character(iter.i)]]);
		sorted.names <- sort(names(nodes.frequency.table), decreasing = TRUE);
		line.plot.data[["degree"]][["x"]] <- nodes.frequency.table[sorted.names];
		line.plot.data[["degree"]][["y"]] <- nodes.frequency.table;

		make.scatter.plot(
			x = as.numeric(names(nodes.frequency.table)),
			y = as.numeric(as.vector(nodes.frequency.table)), 
			correlation.coef = FALSE,
			regression.line = FALSE,
			lowess.line = TRUE,
			file.name = paste(graphs.directory, "nodes_degree__iteration_", iteration, "_", iter.i, ".png", sep = ""),
			height = 5,
			width = 5,
			yaxis.label.offset = -2.91,
			xaxis.label = "Betweeness (fraction)",
			yaxis.label = "Number of nodes",
			resolution = 500
			);
		}

	return (1);
}
