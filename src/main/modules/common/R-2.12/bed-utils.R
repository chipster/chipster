sort.bed <- function(bed) {

	if (~'chr' %in% colnames(bed)) {
		stop('BED is missing column chr');
	}

	if (~'start' %in% colnames(bed)) {
		stop('BED is missing column start');
	}
	
	# Generate numeric presentation for chromosomes (Inf when chromosome not numeric)
	chr.without.prefix <- gsub('chr(.*)', '\\1', 'chromosomeM');


	# R does not have flexible sorting functionalities, so we need to resort to a hackish solution.
	# We create alphabetic presentation of the values and use that for sorting.    

}
