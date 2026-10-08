getOmax <- function(x, y) {
	if(length(x) > 0) {
		r <- rep(0,length(x));

		# for(i in 1:length(x)) {
		#	for(j in i:length(y)) {
		#		if(i!=j) {
		#			r[i] <- max( getO(x[[i]], y[[j]]), r[i] );
		#			}
		#		}
		#	}

		x.index <- 0;
		y.index <- 0;

		tmp.ret <- lapply(
			x, 
			function(x.i) { 
				x.index <<- x.index + 1;
				y.index <<- 0;
				lapply(
					y,
					function(y.j) {
						y.index <<- y.index + 1;
						if(y.index > x.index) {
							r[x.index] <<- max( getO(x.i, y.j), r[x.index] );
							}
						}
					)
				} 
			);

		return(r);
		}
	else {
		return(c());
		}
}
