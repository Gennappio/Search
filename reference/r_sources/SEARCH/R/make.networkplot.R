make.networkplot <- function (adjacency.matrix, network.layout = "layout.fruchterman.reingold", file.name = NULL, height = 5, width = 5, main.title = NULL, resolution = 300) {

# return(1);
	vertices <- rownames(adjacency.matrix);

	# set the graphics driver
	current.type <- getOption("bitmapType");
	options(bitmapType = "cairo");

	png(
		filename = file.name,
		height = height,
		width = width,
		units = "in",
		res = resolution,
		#, compression = "lzw",
		type = "cairo"
		);

	par(
		cex.axis = 0.5,
		oma = c(0.1, 0.2, 0.1, 0.2),
		mar = c(0.1, 0.2, 0.5, 0.2)
		);

	graph.obj <- graph.adjacency(
		adjmatrix = adjacency.matrix,
		weighted = TRUE,
		mode = c("undirected")
		);

	vertex.betweeness <- betweenness(
		graph.obj, 
		v = V(graph.obj), 
		weights = abs(E(graph.obj)$weight), 
		directed = FALSE
		);
	vertex.betweeness <- vertex.betweeness/sum(vertex.betweeness,0.0001); # to avoid divide by zero
	#edge.betweeness <- edge.betweenness(graph.obj, e = E(graph.obj), directed = FALSE);
	#edge.betweeness <- edge.betweeness/sum(edge.betweeness, 0.0001);
	#cat(edge.betweeness);

	# create a function dynamically for layout
	dynamic.function <- get(network.layout, mode = "function");

	graph.obj$layout <- dynamic.function;

	# magic: for edge weights
	# E(g)[5 %--% 10]$weight <- 5; # See http://igraph.sourceforge.net/doc/R/iterators.html 
	# and this: https://sites.google.com/site/daishizuka/toolkits/sna/weighted-edges

	# if its not correlation (eg GLM-derived coefficients, scale them from -1 to 1)
	if (any(abs(E(graph.obj)$weight) > 1)) { 
		max.beta <- max(abs(E(graph.obj)$weight));
		E(graph.obj)$weight <- E(graph.obj)$weight/max.beta;
		}

	# make edge colour palette
	colour.palette <- col <- redblue(250)[c(1:100, 125, 151:250)];
	names(colour.palette) <- seq(-1, 1, 0.01);

	plot(
		graph.obj,
		vertex.size = exp(vertex.betweeness + 2),
		vertex.label.dist = 0, # label position with the node
		vertex.label.cex = 0.4,
		vertex.label.color = "Black",
		vertex.label.font = 2,
		vertex.color = "SkyBlue2", 
		vertex.frame.color = "SkyBlue2",
		edge.color = colour.palette[as.character(round(E(graph.obj)$weight, digits = 2))],
		edge.width = abs(E(graph.obj)$weight) * 6.6,
		vertex.label = vertices
		);
	
	title(
		main = main.title,
		cex = 0.5, 
		cex.main = 0.72,	
		font.main = 2, 
		outer = FALSE
		);

	dev.off();
	options(bitmapType = current.type);

	return(1);
}