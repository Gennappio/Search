nodes.glm <- function(seeds, x, Blist, covariates = NULL, significance.threshold = 0.01, adjust.p = NULL, glm.family = "gaussian", formula.string = NULL, mc.cores = 1) {

	# stop();
	glm.families.P.col <- c(
		"gaussian" = "Pr(>|t|)",
		"Gamma" = "Pr(>|t|)",
		"inverse.gaussian" = "Pr(>|t|)",
		"poisson" = "Pr(>|z|)",
		"quasi" = "Pr(>|t|)",
		"quasipoisson" = "Pr(>|t|)"
		);

	# output list
	glm.res.all <- list();
	glm.res.beta <- list();

	if(glm.family %in% names(glm.families.P.col)) {

		x <- x[rownames(x) %in% Blist, ];

		# create model formula depending upon covariates and formula.string
		if (!is.null(covariates) && !is.null(formula.string)) {
			model.formula <- as.formula(
				paste("seed.exp ~ x.node + ", formula.string)
				);
			all.data <- data.frame(
				"x.node" = rep(NA, ncol(x)), 
				covariates, 
				stringsAsFactors = FALSE
				);
			}
		else {
			model.formula <- as.formula( "seed.exp ~ x.node" );
			all.data <- data.frame(
				"x.node" = rep(NA, ncol(x)), 
				stringsAsFactors = FALSE
				);
			}

		# store rownames for quick lookup
		x.rownames <- rownames(x);

		# go through all seeds one by one
		for (seed in seeds) {

			seed.exp <- x[seed, ];

			# fit glm, and retrieve beta and P value
			glm.res <- mclapply(
				x.rownames,
				FUN = function(x.rowname, x.matrix) {
					x.node <- x.matrix[x.rowname, ];
					all.data$x.node <- x.node;

					# fitting GLM
					glm.fit.obj <- glm(
						model.formula, 
						family = glm.family, 
						data = all.data
						)

					# as opposed to using summary(), do the following to cater for missing rows
					# from coef when due to multicollinearity (singularity) certain interactions
					# are dropped from summary() output. below adds NA.
					glm.fit.coef <- coef(summary(glm.fit.obj));
					glm.fit.tmp <- cbind(
						"Estimate" = coef(glm.fit.obj),
						"P" = rep(NA, length(coef(glm.fit.obj))) 
						);
					glm.fit.tmp[rownames(glm.fit.coef), "P"] <- glm.fit.coef[, glm.families.P.col[glm.family]];

					# extract estimates
					glm.beta <- glm.fit.tmp["x.node", "Estimate"];
					glm.P <- glm.fit.tmp["x.node", "P"];

					# rename row and colnames for unique row and colnames
					rownames(glm.fit.tmp) <- paste(seed, rownames(glm.fit.tmp), sep = ".");
					colnames(glm.fit.tmp) <- paste(x.rowname, colnames(glm.fit.tmp), sep = ".");

					return(
						list(
							"beta" = glm.beta,
							"P" = glm.P,
							"coef.table" = glm.fit.tmp
							)
						)
					},
				x.matrix = x,
				mc.cores = mc.cores
				);

			cat(".");
			# recover some memory
			for(i in 1:5) { gc(); }

			names(glm.res) <- x.rownames;
			glm.res.table <- do.call(rbind, glm.res);
#print(str(glm.res));
#print(glm.res);
#print(glm.res.table);
#stop("\n");

			# correct P values if specified
			if (!is.null(adjust.p)) {
				glm.res.table[, "P"] <- p.adjust(glm.res.table[, "P"], method = adjust.p);
				}
			glm.res.filtered <- glm.res.table[which(glm.res.table[, "P"] < significance.threshold), "beta"];
			glm.res.filtered <- unlist(glm.res.filtered);
			glm.res.beta[[seed]] <- glm.res.filtered;

			# store full coefficients table for interpretation. This include all models,
			# not just the ones which satisfy beta and P thresholds
			glm.res.all[[seed]] <- NULL;
			tmp.obj <- mclapply(
				glm.res, 
				FUN = function(list.element) { list.element$coef.table },
				mc.cores = mc.cores
				);
			glm.res.all[[seed]] <- do.call(cbind, tmp.obj);
			rm(tmp.obj);
			}

#cat("\nGET CORRELATED - GLM:\n");
#print(length(glm.res.all));
#for(i in 1:length(glm.res.all)) {
#	print(dim(glm.res.all[[i]]));
#	}
#print(str(glm.res.filtered));
#print(glm.res.filtered);
#cat("\n");

		# flatten coefficients table
		glm.res.all <- do.call(rbind, glm.res.all);

		return( 
			list(
				"CL" = glm.res.beta,
				"other.data" = glm.res.all
				) 
			);
		}

	stop("glm.family not supported!");
	}
