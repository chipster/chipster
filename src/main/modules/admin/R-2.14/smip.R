#
# SMart Install Packages (SMIP)
#
# User friendly interface into install.packages(...) and biocLite(...)
#

# Usage examples:
#
#smart.install.bioconductor.group("lite")
#smart.install.bioconductor.repo(repo.index = 4)
#smart.install.packages(package="Matrix")
#smart.install.packages(package="pökälestatistics") # should fail
#smart.install.packages(bioconductor.package="Matrix")
#smart.install.packages(bioconductor.package="biopökäle") # should fail
#smart.install.packages(package="Matrix", bioconductor.package="biopökäle") # should fail
#smart.install.packages(url.package="http://cran.r-project.org/src/contrib/00Archive/hdrcde/hdrcde_2.15.tar.gz")


smart.install.bioconductor.group <- function(groupName = "lite", mirror=NA) {
	for (package in biocinstallPkgGroups(groupName)) {
		smart.install.packages(bioconductor.package = package, mirror = mirror)
	}
}

smart.install.bioconductor.repo <- function(repo.index, mirror=NA) {
	
	# Select the given repository
	orig.repos <- setRepositories(ind=c(repo.index))
	
	for (package in available.packages()) {
		smart.install.packages(bioconductor.package = package[1], mirror = mirror)
	}
	
	# Restore original repositories
	setRepositories(orig.repos)		
}

smart.install.packages <- function(package=NA, bioconductor.package=NA, url.package=NA, mirror=NA) {
	
	# Check parameters
	package.defs <- c(package, bioconductor.package, url.package)
	if (sum(!is.na(package.defs)) != 1) {
		stop("Must use exactly one of the alternative ways to specify the package to install")
	}

	# Check that not already installed		
	if (is.na(url.package)) {
		package.name <- package.defs[!is.na(package.defs)]
	} else {
		package.name <- gsub(".*/(.*)_.*", "\\1", url.package)
	}
	
	if (!is.installed(package.name)) {
		cat(paste("Will now install", package.name, "\n"))
	} else {
		cat(paste("Already installed", package.name, "\n"))
		return(TRUE)
	}
	
	if (!is.na(package)) {
		
		repos = ifelse(is.na(mirror), getOption("repos"), mirror)
		install.packages(packages=c(package), repos)
	
	} else if (!is.na(bioconductor.package)) {
		
		source("http://www.bioconductor.org/biocLite.R")
		if (!is.na(mirror)) {
			options("BioC_mirror" = c("Mirror"=mirror))
		}
		biocLite(bioconductor.package)
		
	} else if (!is.na(url.package)) {
		
		# Download URL to temp file, install and remove
		tempfile.path <- tempfile("package", fileext=".tar.gz")
		download.file(url=url.package, destfile=tempfile.path)
		install.packages(packages=c(tempfile.path))
		unlink(tempfile.path)
		
		
	} else {
		stop("Must specify something to install");
	}

	# Install was successfull
	return(TRUE)
}

is.installed <- function(package) {
	
	return(suppressWarnings(require(package, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE))) 		
}
