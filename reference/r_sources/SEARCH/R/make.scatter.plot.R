make.scatter.plot <- function (x = NULL, y = NULL, correlation.method = "spearman", correlation.coef = TRUE, regression.line = TRUE, lowess.line = TRUE, xlimit = NULL, x.abline.at = NULL, y.abline.at = NULL, cex.axis = 1.2, cex.lab = 1.4, main.title = "", xaxis.label = expression(bold(paste('log'[2], ' mRNA abundance'))), yaxis.label = expression(bold(paste('log'[2], ' mRNA abundance'))), yaxis.label.offset = -1.4, margins = c(2.7, 3.5, 0.5, 0.2), file.name = NULL, height = 5, width = 5, resolution = 300) {

	# sanity checks
	mar.bottom <- margins[1];
	mar.left <- margins[2];
	mar.top <- ifelse (main.title != "" && !is.null(main.title), 1.2, margins[3]);
	mar.right <- margins[4];
	mgp.lab <- 1.66;
	if (length(x) < 5) {
		cat("\nLess than 5 data points. Skipping correlation, regression and lowess stats on scatter plot");
		correlation.coef <- FALSE; regression.line <- FALSE; lowess.line <- FALSE;
		}

	# set the graphics driver
	current.type <- getOption("bitmapType");
	options(bitmapType = "cairo");

	png(
		filename = file.name,
		height = height,
		width = width,
		units = "in",
		res = resolution
		#, compression = "lzw"
		);

	par(
		tcl = -0.2,
		las = 1,
		font.lab = 2,
		cex.lab = cex.lab,
		font.axis = 2,
		cex.axis = cex.axis,
		mar = c(mar.bottom, mar.left, mar.top, mar.right),
		mgp = c(mgp.lab, 0.25, 0),
		oma = c(0, 0, 0, 0),
		family = "sans"
		);

	plot(
		x, 
		y,
		xaxs = "r",
		xlab = xaxis.label,
		ylab = "",
		pch = 20,
		lwd = 2
		);

	# yaxis label
	mtext(
		las = 0,
		padj = yaxis.label.offset, # distance from the y-axis line
		text = yaxis.label,
		side = 2,
		font = 2,
		cex = cex.lab,
		font.axis = 2
		);

	# main title
	title(main.title, cex.main = cex.lab, font.main = 2);

	# regression line
	if (regression.line) {
		abline(lm(y ~ x), col = "red") # regression line (y~x) 
		}

	# Lowess's fit
	if (lowess.line) {
		lines(lowess(x, y), col = "blue") # lowess line (x,y)
		}

	# print correlation rho
	# par(mar = c(0, 0, 0, 0), font = 2);
	if (correlation.coef) {

		# correlation estimate
		cor.rho <- cor.test(x, y, method = correlation.method)$estimate;

		corner.label(
			expression(rho),
			x = -1, y = 1, 
			figcorner = FALSE
			);
		corner.label(
			paste("=", round(cor.rho, digits = 2)),
			x = -1, y = 1, 
			xoff = 0.086,
			figcorner = FALSE
			);
		}

	if (regression.line && lowess.line) {
		par(mar = c(0, 0, 0, 0), font = 2);
		legend(
			x = "topleft", 
			legend = c("Regression", "Lowess"),
			col = c("red", "blue"),
			lwd = 2,
			bty = "n",
			cex = 0.9,
			adj = c(0, 0.203),
			inset = c(0, 0.05)
			);
		}

	# make background grid if requested
	if (length(x.abline.at) > 0) { abline(v = x.abline.at, col = "grey", lty = "dashed") };
	if (length(y.abline.at) > 0) { abline(h = y.abline.at, col = "grey", lty = "dashed") };

	dev.off();
	options(bitmapType = current.type);

	# return true
	return ( 1 );

	}
