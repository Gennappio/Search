# removeNode <- function (..., list = character(), pos = -1, envir = as.environment(pos), inherits = FALSE) {
#	dots <- match.call(expand.dots = FALSE)$...
#	if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || is.character(x)))) 
#		stop("... must contain names or character strings")

#	names <- sapply(dots, as.character)
#	if (length(names) == 0L) 
#		names <- character()

#	list <- .Primitive("c")(list, names)
#	.Internal(remove(list, envir, inherits))
#}

removeNode <- function(to.remove, x) {

	if (length(to.remove) > 0) {
		to.remove.index <- which(x %in% to.remove);
		if (length(to.remove.index) > 0) {
			x <- x[-to.remove.index];
			}
		}

	return(x);
}
