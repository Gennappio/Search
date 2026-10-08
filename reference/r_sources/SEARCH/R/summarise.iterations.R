summarise.iterations <- function(output.directory = ".", max.iter = 1, n.boot = 1, output.directory.NULL = NULL, n.boot.NULL = NULL, batches.NULL = NULL) {

	# initialise variables
	features.all <- vector();
	conn.score.mean <- list();
	conn.score.mean.NULL <- list();

	# compile data for all iterations
	for (iter in as.character(1:max.iter)) {

		conn.score.all <- list();
		conn.score.all.NULL <- list();

		# go through [REAL] runs and compile connectivity scores
		for (i in as.character(1:n.boot)) {

			# read connectivity score and take mean across iterations
			conn.file <- paste(
				output.directory, 
				"output/connectivity__iteration_", i, "_", iter, ".txt.gz", 
				sep = ""
				);

			if (file.exists(conn.file)) {
				conn.score.all[[i]] <- as.matrix(
					read.table(
						file = conn.file,
						row.names = 1,
						header = TRUE,
						sep = "\t"
						)
					);

				features.all <- union(features.all, rownames(conn.score.all[[i]]));
				}
			else {
				cat("\nskipping file (does not exists): ", conn.file);
				}
			}

		# initialise and store all iterations systematically
		conn.score.mean[[iter]] <- matrix(
			data = NA, 
			nrow = length(features.all),
			ncol = length(conn.score.all),
			dimnames = list(
				features.all,
				1:length(conn.score.all)
				)
			);

		for (i in colnames(conn.score.mean[[iter]])) {
			conn.score.mean[[iter]][rownames(conn.score.all[[i]]), i] <- conn.score.all[[i]][, 1];
			}

		# add mean column
		conn.score.mean[[iter]] <- cbind(
			conn.score.mean[[iter]], 
			"mean" = rowMeans(conn.score.mean[[iter]], na.rm = TRUE)
			);

		# order by mean value
		conn.score.mean[[iter]] <- conn.score.mean[[iter]][
			order(conn.score.mean[[iter]][, "mean"], decreasing = TRUE), 
			];

		# go through [NULL] runs and compile connectivity scores
		# process data from NULL distributions to estimate a P value and add that the P value column
		# to the conn.score.mean
		conn.score.all.NULL <- list();
		for (batch.i in 1:batches.NULL) {
			for (i in as.character(1:n.boot.NULL)) {

				# read connectivity score and take mean across iterations
				conn.file <- paste(
					paste(output.directory.NULL, "batch_", batch.i, "_0_FALSE_0_0/", sep = ""), 
					"output/connectivity__iteration_", i, "_", iter, ".txt.gz", 
					sep = ""
					);

				if (file.exists(conn.file)) {
					i.effective <- as.numeric(i) + (n.boot.NULL * (batch.i - 1));
					#cat("\ni NULL index: ", i.effective);
					conn.score.all.NULL[[as.character(i.effective)]] <- as.matrix(
						read.table(
							file = conn.file,
							row.names = 1,
							header = TRUE,
							sep = "\t"
							)
						);

					#print(head(conn.score.all.NULL[[i.effective]]));
					}
				else {
					cat("\nskipping file (does not exists): ", conn.file);
					}
				}
			}

		# initialise p-values column
		conn.score.mean[[iter]] <- cbind(
			conn.score.mean[[iter]],
			"P" = rep(NA, nrow(conn.score.mean[[iter]]))
			);

		# if data on NULL distributions available
		if (length(conn.score.all.NULL) > 0) {

			# initialise and store all iterations systematically (NULL distributions)
			conn.score.mean.NULL[[iter]] <- matrix(
				data = NA, 
				nrow = length(features.all),
				ncol = length(conn.score.all.NULL),
				dimnames = list(
					features.all,
					1:length(conn.score.all.NULL)
					)
				);
			# sort as the non-NULL data
			conn.score.mean.NULL[[iter]] <- conn.score.mean.NULL[[iter]][rownames(conn.score.mean[[iter]]), ];

			for (i in colnames(conn.score.mean.NULL[[iter]])) {
				conn.score.mean.NULL[[iter]][rownames(conn.score.all.NULL[[i]]), i] <- conn.score.all.NULL[[i]][, 1];
				}

			# estimate P value for each gene/feature
			for (feature.i in rownames(conn.score.mean[[iter]])) {
				conn.score.mean[[iter]][feature.i, "P"] <- (
					length(which(conn.score.mean.NULL[[iter]][feature.i, ] >= conn.score.mean[[iter]][feature.i, "mean"])) + 1
					)/(ncol(conn.score.mean.NULL[[iter]]) + 1);
				}

			# order by P-value
			# conn.score.mean[[iter]] <- conn.score.mean[[iter]][
			#	order(conn.score.mean[[iter]][, "P"], decreasing = FALSE), 
			#	];
			}

		} # end of iter loop

	return(list("real" = conn.score.mean, "NULL" = conn.score.mean.NULL));
	}
