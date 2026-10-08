removeByName <- function(x, l) {
	# i.e. keep everything whos name isn't in 'x'
	l[!(names(l) %in% x)] 
}
