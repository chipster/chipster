unzipIfGZipFile <- function(file.name){
	
	# get file type with the unix file command
	file.type = system(paste("file -Lb", file.name), intern=TRUE)

	# if gzip, unzip it, TODO fix that startsWITH
	if (startsWith(filetype, "gzip")) {
		zipfile.name <- paste(file.name, ".gz", sep="")
		system(paste("mv", filename, zipfile.name, "; gzip -d", zipfile.name))
	}
}