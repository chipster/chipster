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

smart.install.scavenge.web.packages <- function(url, chiptype="all", update=0, list.only=0) {
	# Name of the static pages
	if(length(grep("latest", url)) > 0) {
		web.page <- paste(readLines(url), collapse="")
		latest.build <- gsub("\'.*$", "", gsub("^.*URL=\\.\\.", "", web.page))
		url <- gsub("\\/latest\\/.*$", latest.build, url);
	}

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

	if(chiptype != "all") {
		links <- links[grep(chiptype, links)];
	}

	# if list.only == 1, only the list of available packages will returned
	if(list.only == 1) {
		return(links)
	}

	# Install each linked package
	for (link in links) {
		smart.install.packages(url.package=link, update=update)
	}
}

smart.install.bioconductor.repo <- function(repo.index, mirror=NA, update=0) {
	
	# check what version of Bionconductor is being used
	bioc.intaller.loaded <- try(library(BiocInstaller))
	if(class(bioc.intaller.loaded) == "try-error") {
		smart.install.packages(bioconductor.package="BiocInstaller", mirror=mirror, update=update)
	}
	library(BiocInstaller)
	current.bioc.version <- biocVersion()

	# Install available annotation packages
	current.bioc.url <- paste("http://www.bioconductor.org/packages/", current.bioc.version, "/data/annotation", sep="");
	for (package in available.packages(contrib.url(current.bioc.url))[,"Package"]) {
		smart.install.packages(bioconductor.package = package, mirror = mirror, update=update)
	}

	# Select the given repository
	#orig.repos <- setRepositories(ind=c(repo.index))
	#
	#for (package in available.packages()[,"Package"]) {
	#	smart.install.packages(bioconductor.package = package, mirror = mirror, update=update)
	#}
	#
	# Restore original repositories
	#if(length(orig.repos) == 1) {
	#	if(orig.repos$repos["CRAN"] == "@CRAN@") {
	#		#orig.repos$repos["CRAN"] = "@CRAN@";
	#		options(repos = orig.repos);
	#	} else {
	#		setRepositories(orig.repos)
	#	}
	#} else {
	#	setRepositories(orig.repos)
	#}
}

smart.install.packages <- function(package=NA, bioconductor.package=NA, url.package=NA, mirror=NA, update=0) {
	
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

	if(update==0) {
		if (!is.installed(package.name)) {
			cat(paste("Will now install", package.name, "\n"))
		} else {
			cat(paste("Already installed", package.name, "\n"))
			return(invisible(TRUE))
		}
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
		a <- try(download.file(url=url.package, destfile=tempfile.path))
		if(class(a) != "try-error") {
			install.packages(pkgs=c(tempfile.path), repos=NULL)
			unlink(tempfile.path)
		} else {
			warning(paste("package", url.package, "is not valid a web-page", sep=" "))
		}
		
	} else {
		stop("Must specify something to install");
		return(invisible(FALSE))
	}

	# Install was successfull
	return(invisible(TRUE))
}

is.installed <- function(package) {
	if(package %in% rownames(installed.packages()) == FALSE) {
		sink("/dev/null") # the only way to get rid of all output (some packages don't behave)
		is.installed <- suppressPackageStartupMessages(suppressWarnings(suppressMessages(require(package, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE))))
		sink()
		return(is.installed)
	} else {
		return(package %in% rownames(installed.packages()))
	}
}

check.affy.customnames <- function(script.basename, db.custom.packages) {

	# Find R scripts used to normalize affy chips
	r.script.basename <- gsub("/admin/", "/microarray/", script.basename);
	all.packages <- rownames(installed.packages())
	all.r.scripts <- list.files(path=r.script.basename);
	affy.r.script <- all.r.scripts[grep("norm-affy", all.r.scripts)];

	# Install stringdist package. Needed for finding partial matches
	smart.install.packages(package="stringdist", mirror=repo.cran)
	library("stringdist");

	# Polish up custom-CDF names
	if(length(db.custom.packages) > 0) {
		db.custom.packages <- gsub(".*/(.*)_.*", "\\1", db.custom.packages)
	} else {
		stop("Must specify a list of custom-CDF packages");
	}

	# For each affy-norm script, find code defining which custom-CDF packages are supported
	supported.custom.packages <- NULL;
	for(i in 1:length(affy.r.script)) {
		r.script <- scan(file=file.path(r.script.basename, affy.r.script[i]), what="", sep="\n")

		#Find instances in SADL descriptions
		if(length(grep("PARAMETER custom\\.chiptype", r.script) > 0)) {
			sadl.row <- r.script[grep("PARAMETER custom\\.chiptype", r.script)];
			sadl.row <- gsub("^.*\\[", "", sadl.row);
			sadl.row <- gsub("\\].*$", "", sadl.row);
			packages.in.sadl.row <- unlist(strsplit(sadl.row, "\\s*,\\s*"))
			for(j in 1:length(packages.in.sadl.row)) {
				custom.package <- unlist(strsplit(packages.in.sadl.row[j], "\\s*:\\s*"))[2];
				custom.package <- gsub("\\(.+\\)", "cdf", custom.package);
				supported.custom.packages <- c(supported.custom.packages, custom.package);
			}
		}

		#Find other instances where parameter custom_cdf has been used
		if(length(grep("custom_cdf\\s*<-|custom_cdf\\s*=", r.script) > 0)) {
			rscript.row <- r.script[grep("custom_cdf\\s*<-|custom_cdf\\s*=", r.script)];
			rscript.row <- gsub("^.*<-|^.*=", "", rscript.row);
			rscript.row <- gsub("\\s+|\"", "", rscript.row);
			supported.custom.packages <- c(supported.custom.packages, rscript.row);
		}
	}

	# Check if the package exists
	for(j in 1:length(supported.custom.packages)) {
		if(!(supported.custom.packages[j] %in% db.custom.packages)) {
			cat(paste("Package", supported.custom.packages[j], "in", affy.r.script[i], "not found\n"));
			partial.match <- db.custom.packages[ain(db.custom.packages, supported.custom.packages[j], maxDist=3)];
			if(length(partial.match) > 0) {
				for(k in 1:length(partial.match)) {
					if(partial.match[k] %in% rownames(installed.packages()) == TRUE) {
						cat(paste("\tConsider using", partial.match[k], "\n"));
					} else {
						cat(paste("\tConsider installing and using", partial.match[k], "\n"));
					}
				}
			} else {
				cat(paste("\tPackage", supported.custom.packages[j], "in", affy.r.script[i], "has not matches in current custom-CDF database\n"));
			}
		}
	}
}

