
mergeFiles <- function(files) {
	print(paste0("Total files found: ", length(files)))
	df <- NULL
	for (i in seq(files)) {
	    print(paste0("Reading file: ", i, "/", length(files)))
	    f <- files[i]

	    if (is.null(df)) {
	        df <- read.table(f, sep ="\t", header=TRUE)
	    } else {
	        temp <- read.table(f, sep ='\t', header=TRUE)
	        df <- rbind(df, temp)
	    }
	}

	print(paste0("Total of data imported: ", dim(df)[1]))
	return(df)
}