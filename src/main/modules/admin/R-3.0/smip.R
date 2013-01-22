#
# SMart Install Packages (SMIP)
#
# User friendly interface into install.packages(...) and biocLite(...)
#
# All functions check existence of packages before installation and skip installation
# if the package already exists. For group of packages, each one is checked
# individually and only missing ones are installed.
#
# Function also offer more automation compared to original ones, allowing
# installation of whole repositories, scavenging web pages for packages etc.

# Usage examples:
#
#smart.install.bioconductor.repo(repo.index = 4)
#smart.install.packages(package="Matrix")
#smart.install.packages(package="Matrix", mirror="http://ftp.sunet.se/pub/lang/CRAN")
#smart.install.packages(package="foobarstatistics") # should fail
#smart.install.packages(bioconductor.package="Matrix")
#smart.install.packages(bioconductor.package="biofoobar") # should fail
#smart.install.packages(package="Matrix", bioconductor.package="biofoobar") # should fail
#smart.install.packages(url.package="http://www.math.utu.fi/projects/software/bio/ROTS_1.1.1.tar.gz")
#smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/entrezg.asp")

smart.install.scavenge.web.packages <- function(url) {

	# Parse URL
	url.root <- substring(url, 1, gregexpr("/[^/]*$", url)[[1]][1])
	
	# Read page that contains links
	web.page <- paste(readLines(url), collapse="")
	
	# Parse links
	links.raw <- regmatches(web.page, gregexpr("href=\"([^\"]*)\"", web.page))[[1]]
	links.stripped <- substr(links.raw, nchar("href=\"")+1, nchar(links.raw)-1)
	links <- links.stripped[grep(".tar.gz", links.stripped)]
	
	# Convert relative URL's to absolute 
	ind.relative.url <- grep("http://", links, invert=TRUE)
	links[ind.relative.url] <- paste(url.root, links[ind.relative.url], sep="")
	
	# Install each linked package
	for (link in links) {
		smart.install.packages(url.package=link)
	}
}

smart.install.bioconductor.repo <- function(repo.index, mirror=NA) {
	
	# Select the given repository
	orig.repos <- setRepositories(ind=c(repo.index))
	
	for (package in available.packages()[,"Package"]) {
		smart.install.packages(bioconductor.package = package, mirror = mirror)
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
		return(invisible(TRUE))
	}
	
	if (!is.na(package)) {
		
		repos = ifelse(is.na(mirror), getOption("repos"), mirror)
		install.packages(pkgs=c(package), repos=repos)
	
	} else if (!is.na(bioconductor.package)) {
		
		source("http://www.bioconductor.org/biocLite.R")
		if (!is.na(mirror)) {
			options("BioC_mirror" = c("Mirror"=mirror))
		}
		biocLite(bioconductor.package, suppressUpdates=TRUE)
		
	} else if (!is.na(url.package)) {
		
		# Download URL to temp file, install and remove
		tempfile.path <- tempfile("package", fileext=".tar.gz")
		download.file(url=url.package, destfile=tempfile.path)
		install.packages(pkgs=c(tempfile.path), repos=NULL)
		unlink(tempfile.path)
		
		
	} else {
		stop("Must specify something to install");
		return(invisible(FALSE))
	}

	# Install was successfull
	return(invisible(TRUE))
}

is.installed <- function(package) {
	sink("/dev/null") # the only way to get rid of all output (some packages don't behave)
	is.installed <- suppressPackageStartupMessages(suppressWarnings(suppressMessages(require(package, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE))))
	sink()
	return(is.installed) 		
}
