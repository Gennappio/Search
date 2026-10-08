SEARCH <- function(x, Alist, Blist, covariates = NULL, scale = TRUE, D.lower = 0, D.higher = 0, O.lower = 0, keep.Alist = TRUE, max.iter = 1, sd.threshold = 0, verbose = TRUE, nodes.association.function = "nodes.correlation", ...) {

	# initialise variables
	D <- list();
	O <- matrix();
	Omax <- list();
	correlationLists <- list();
	removed <- vector();
	added <- c();
	iter <- 1;
	network.growth <- list();
	network.connectivity <- list();
	network.correlation <- list();
	network.overlap <- list();
	other.data <- list();

	# Make sure Blist starts off with the Alist in it
	Blist <- union(Alist, Blist);

	# pre-filter on standard deviation
	# remove any probesets with an sd across samples less than the threshold
	if(sd.threshold > 0 ) {
		if(verbose) {
			cat("Prefiltering on standard deviation ... ");
			}
		sds <- apply(x[rownames(x) %in% Blist,], 1, sd);
		names(sds) <- Blist;
		to.remove <- names(sds[sds < sd.threshold]);
		Alist <- removeNode(to.remove, Alist);
		Blist <- removeNode(to.remove, Blist);
		if(verbose){
			cat("done.\n");
			cat("Removed ", length(to.remove), " nodes because their standard deviation is below the threshold.\n");
			}
		}

	# filtered seeds and universe list
	seed.nodes.filtered <- Alist;
	universe.nodes.filtered <- Blist;

	# iterate over desired number of segments away from the seed nodes
	while(iter <= max.iter) {

		if(verbose) { cat("\nmax.iter =", iter, "/", max.iter, "\n"); }

		n.remO <- 0;

		# check if the Alist is empty. If it is stop.
		if(length(Alist) == 0) { 
			cat("Stopping, Alist is empty.\n");
			correlationLists <- list();
			break;
			}

		# Use lapply to process the Alist one by one.
		# For each member, get the Blist members it is correlated with
		# the result is a list of lists (see above).

		if(verbose) { cat("\tGetting correlated sets ... "); }

		# Only calculate all correlations if we need to. Otherwise simply do the new stuff
		# this iteration and append the results to the existing list.

		# create a function dynamically
		dynamic.function <- get(nodes.association.function, mode = "function");
		

		if(iter == 1) {

			association.data <- dynamic.function(
				seeds = Alist, 
				x = x,
				Blist = Blist,
				covariates = covariates,
				...
				);

			correlationLists <- association.data[["CL"]];
			}
		else {

			# just do the new stuff *SH: some description needed here*
			if(length(added) > 0) {

				association.data <- dynamic.function(
					seeds = added, 
					x = x,
					Blist = Blist,
					covariates = covariates,
					...
					);

				newCL <- association.data[["CL"]];

				# print(names(newCL));
				# newCL contains stuff added because D was high enough. 
				# Need to check these to see if O is OK...

				# GO THROUGH THESE AND REMOVE BASED ON O
				if(verbose) { cat("\n\tCalculating max(O) ..."); }

				# Alist is not immune to removal (to identify a compact signature)
				if(!keep.Alist) {
					cat("\n\t\t\t\tseeds are NOT immune to removal\n");
					Omax <- getOmax(x = c(newCL, correlationLists), y = c(newCL, correlationLists));
					names(Omax) <- c(added, names(correlationLists));
					}
				else {
					cat("\n\t\t\t\tseeds are immune to removal\n");
					Omax <- getOmax(newCL, c(newCL, correlationLists));
					names(Omax) <- added;
					}

				if(verbose) { cat("\tdone.\n"); }

				to.remove <- Omax;
				to.remove <- names(to.remove[to.remove < O.lower]);
				Alist <- removeNode(to.remove, Alist);

				# note, we take from the Blist too - they're gone permanently
				# *SH: some description needed here*
				Blist <- removeNode(to.remove, Blist);

				newCL <- removeByName(to.remove, newCL);

				# how many did we get rid of
				n.remO <- length(to.remove);

				if(verbose) {
					if(length(to.remove) >  0) {
						cat("\tRemoved the following nodes based on Omax:");
						print(to.remove);
						}
					else {
						cat("\tRemoved nothing based on Omax");
						}
					}
				correlationLists <- c(newCL, correlationLists);
				}
			}

		if(verbose) { cat("\tdone.\n"); }

		# We're going to count how many times each of the probesets 
		# we're considering appears in a correlationList
		# Set up a counter for each probeset in the Blist

		if(verbose) { cat("\tCalculating D ...") }

		counts <- rep(0, length(Blist));
		names(counts) <- Blist;

		# Create a list of lists. 
		# The outer list is the length of the Alist.
		# Each member is the length of the Blist (and in the same order),
		# and contains T if that probeset is in the Alist set, F if not.
		# Then run along that list of lists and add 1 to count every time we see a T.

		toCount <- lapply(
			correlationLists,
			function(tmp.x,all.names) {
				tmp.x <- names(tmp.x)
				all.names %in% tmp.x 
				},
			all.names = Blist
			);

		# we now have a list, one entry for each Alist member. Each of these is a T/F list, one entry
		# for each Blist member. T when the Blist member is in the corresponding Alist member's correlation
		# list.
		# We want to count, for each Blist member, how many times that happens

		for(i in 1:length(toCount)) {
			tmp.x <- toCount[[i]];
			# i.e. for each member of count, where tmp.x is true, add 1
			counts[tmp.x] <- counts[tmp.x] + 1;
			}

		# adjust for self-interactions
		counts[names(toCount)] <- counts[names(toCount)] - 1;

		# we get D by scaling these counts by the number of possible intersections we could have
		if(scale) {
			# adjust for self-interactions
			counts[names(toCount)] <- counts[names(toCount)]/(length(Alist)-1);
			# scale for non-self interactions
			counts[setdiff(Blist, names(toCount))] <- counts[setdiff(Blist, names(toCount))]/length(Alist);
			}
		D <- counts;

		if(verbose) { cat("\tdone.\n"); }

		# So lets remove things from the Alist according to our remove rules:

		# remove nodes based on D
		# remove anything that has a value of D lower than the threshold 'D.lower'
		to.remove <- D[names(D) %in% Alist];
		to.remove <- names(to.remove[to.remove < D.lower]);

		Alist <- removeNode(to.remove, Alist);
		Blist <- removeNode(to.remove, Blist);

		# note, we take from the Blist too - they're gone permanently
		correlationLists <- removeByName(to.remove, correlationLists);

		# don't forget D and Omax
		D <- removeByName(to.remove, D);
		Omax <- removeByName(to.remove, Omax);

		network.connectivity[[as.character(iter)]] <- D;

		removed <- c(to.remove, removed);

		if(verbose) {
			if(length(to.remove) > 0) {
				cat("\tRemoved the following nodes based on D:\n");
				print(to.remove);
				}
			else {
				cat("\tRemoved nothing based on D.\n");
				}
			}

		# how many did we get rid of
		n.remD <- length(to.remove);

		# And now we must add probesets to the Alist according to our add rules:
		# ADD TO ALIST BASED ON D
		# and add to Alist, anything where D higher than the threshold 'D.higher'

		not.already <- !(names(D) %in% Alist);
		to.add <- D[not.already];
		to.add <- names(to.add[to.add > D.higher]);

		Alist <- c(Alist, to.add);

		if(verbose) {
			if(length(to.add) > 0) {
				cat("\tAdded the following nodes to the Alist based on D, max(D) =", max(D), "\n");
				print(to.add);
				}
			else {
				cat("\tAdded nothing based on D, max(D) =", max(D), "\n");
				print(to.add);
				}
			}

		# how many did we add *SH: need some description*
		n.addD <- length(to.add);

		# and remember what we added - so we only calculate correlations for the new 
		# stuff in the next iteration
		added <- to.add;
		network.growth[[as.character(iter)]] <- to.add;

		other.data[[as.character(iter)]] <- association.data[["other.data"]];

		# stop if nothing has changed
		# possible return from the function
		if( (n.addD + n.remD + n.remO) == 0 || 
			length(Blist) == 0 ||
			length(correlationLists) == 0
			) {
			if(verbose) {
				cat("\nStopping. Nothing changed in this iteration OR Blist empty, length Blist =", length(Blist)); 
				}
			break;
			}

		cat("\nConverting to Correlation and O matrices...");
		cat("\nsize of CorrelationLists = ", length(correlationLists));
		matrix.format <- get.matrix.format(x = correlationLists, all.features = rownames(x));
		network.correlation[[as.character(iter)]] <- matrix.format[["subnetworks.correlation"]];
		network.overlap[[as.character(iter)]] <- matrix.format[["O"]];

		# increment loop
		iter <- iter + 1;

		# recover some memory
		for(i in 1:5) { gc(); }

		}

	cat("\nDone.\n");

	return(
		list(
			# Alist = Alist, 
			# Blist = Blist, 
			removed = removed, 
			correlationLists = correlationLists, 
			correlationMatrix = network.correlation,
			D = network.connectivity, 
			O = network.overlap,
			network.growth = network.growth,
			other.data = other.data
			)
		);
}
