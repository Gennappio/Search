get.program.defaults <- function(entry.secret = "hypoxia_angiogenesis_sub_networks_all.txt") {

	# make a list of potential locations for the datasets file
	program.data.dirs <- paste(.libPaths(), "/SEARCH/programdata/", sep = "");

	# then search all locations
	file.checks <- file.exists( paste(program.data.dirs, "/", entry.secret, sep = "") );

	# check to see if the file was actually found
	if (any(file.checks)) {
		program.data.dir <- program.data.dirs[ order(file.checks, decreasing = TRUE)[1] ];
		}
	else {
		stop( paste("Unable to find file: ", entry.secret, sep = "") );
		}
	
	return (
		list(
			"program.data.dir" = program.data.dir,
			"subnets.file" = paste(program.data.dir, "/hypoxia_angiogenesis_sub_networks_all.txt", sep = ""),
			"test.data.dir" = paste(program.data.dir, "/testdata/", sep="")
			)
		);
	}
