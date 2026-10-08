getO <- function(x, y) {
	xn <- names(x);
	yn <- names(y);
	length(intersect(xn, yn)) / length(union(xn, yn));
}
