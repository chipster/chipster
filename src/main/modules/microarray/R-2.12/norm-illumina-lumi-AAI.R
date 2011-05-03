# TOOL norm-illumina-lumi-AAI.R: "Illumina - lumi pipeline AAI" (Illumina normalization using BeadSummaryData files that reports the probe identifiers as Array_Address_Id, for example 4900685, instead of Probe_Id, for example ILMN_2607609, and using lumi methodology. TO USE THIS TOOL, YOU NEED TO IMPORT THE BeadSummaryData FILE DIRECTLY, NOT USING THE IMPORT TOOL. Please note that processing the data with Array_Address_Id is slower and the normalization therefore takes a bit longer time than the Illumina - lumi pipeline tool.)
# INPUT chip.tsv: chip.tsv TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER transformation: transformation TYPE [none: none, vst: vst, log2: log2] DEFAULT log2 (How to transform the data)
# PARAMETER background.correction: background.correction TYPE [none: none, bgAdjust.affy: bgAdjust.affy] DEFAULT none (Should background adjustment be applied)
# PARAMETER normalize.chips: normalize.chips TYPE [none: none, rsn: rsn, loess: loess, quantile: quantile, vsn: vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER chiptype: chiptype TYPE [empty: empty, Human: Human, Mouse: Mouse, Rat: Rat] DEFAULT empty ()

# Illumina data preprocessing and normalization for BeadSummaryData that
# uses Array_Address_Id instead of Probe_Id 
#
# MG, 5.11.2010
# NG, 21.12.2010 modified to add gene symbol and gene name to the output

# Loading libraries
library(lumi)
library(annotate)

# Converting to the correct chiptype
if(chiptype=="empty") {
	chiptype<-c("Illumina")
	mapping<-c("Illumina")
}
if(chiptype=="Human") {
	chiptype<-c("lumiHumanAll")
	mapping<-c("lumiHumanIDMapping")
}
if(chiptype=="Mouse") {
	chiptype<-c("lumiMouseAll")
	mapping<-c("lumiMouseIDMapping")
}
if(chiptype=="Rat") {
	chiptype<-c("lumiRatAll")
	mapping<-c("lumiRatIDMapping")
}
chiptype<-paste(chiptype, ".db", sep="")

# Loading data files
# fileName <-dir()
# fileName<-fileName[fileName!="phenodata.tsv"]
x.lumi <- lumiR("chip.tsv", lib.mapping=mapping)

`getChipInfo` <-
		function(x, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), chipVersion=NULL, idMapping=FALSE, returnAllMatches=FALSE, verbose=TRUE) {
	
	## Function to make sure the output mapping matches the inputID.bak
	matchInputID <- function(mapping, inputID.bak) {
		## We assume mapping IDs are a subset of the names of inputID.bak
		nc <- ncol(mapping)
		if (is.null(nc)) {
			len <- length(mapping)
			if (len == length(inputID.bak)) return(mapping)
			mapping.new <- rep(NA, length(inputID.bak))
			names(mapping.new) <- inputID.bak
			mappingID <- names(mapping)
			mapping.new[mappingID] <- mapping
			type <- 'vector'
		} else {
			nr <- nrow(mapping)
			if (nr == length(inputID.bak)) return(mapping)
			mapping.new <- matrix(NA, nrow=length(inputID.bak), ncol=nc)
			rownames(mapping.new) <- inputID.bak
			colnames(mapping.new) <- colnames(mapping)
			mappingID <- rownames(mapping)
			mapping.new[mappingID,] <- mapping[mappingID,]
			type <- 'matrix'
		}
		selID <- inputID.bak[inputID.bak %in% mappingID]
		dupInd <- which(duplicated(selID))
		if (length(dupInd) > 0) {
			dupID <- unique(selID[dupInd])
			for (dupId.i in dupID) {
				ind.i <- which(inputID.bak == dupId.i)
				if (type == 'matrix') {
					mapping.new[ind.i, ] <- mapping[dupId.i, ]
				} else {
					mapping.new[ind.i] <- mapping[dupId.i]					
				}
			}
		}
		return(mapping.new)
	}
	
	
	if (is(x, 'ExpressionSet')) {
		inputID <- featureNames(x)
	} else if (is(x, 'matrix') || is(x, 'data.frame')) {
		inputID <- rownames(x)		
	} else {
		inputID <- as.character(x)
	}
	inputID.bak <- inputID
	inputID <- unique(inputID)
	species <- match.arg(species)
	if (species == 'Unknown' && is.null(lib.mapping)) {
		human.match <- getChipInfo(inputID, species='Human', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		mouse.match <- getChipInfo(inputID, species='Mouse', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		rat.match <- getChipInfo(inputID, species='Rat', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		bestChip <- which.max(c(length(human.match$matchedProbeNumber), length(mouse.match$matchedProbeNumber), length(rat.match$matchedProbeNumber)))
		best.match <- switch(bestChip, human.match, mouse.match, rat.match)
		return(best.match)
	} else {
		if (is.null(lib.mapping)) {
			lib.mapping <- switch(species,
					'Rat'='lumiRatIDMapping',
					'Human'='lumiHumanIDMapping',
					'Mouse'='lumiMouseIDMapping')			
		}		
		if(!require(lib.mapping, character=TRUE)) stop(paste(lib.mapping, 'is required!'))
		# dbconn <- sub("\\.db", "_dbconn", lib.mapping)
		dbconn <- paste(lib.mapping, "_dbconn", sep="")
		conn <- do.call(dbconn, list())
		
		metaInfo <- dbReadTable(conn, 'metadata')
		species <- metaInfo$value[metaInfo[,'name'] == "SPECIES"]
		
		allTableNames <- dbListTables(conn)
		allTableNames <- allTableNames[!(allTableNames %in% c('nuID_MappingInfo', 'metadata'))]
		if (!is.null(chipVersion)) {
			if (chipVersion %in% allTableNames) {
				allTableNames <- chipVersion
			} else {
				warning('The provided chipVersion does not exist in the library!\n')
			}
		}
		lenID <- length(inputID)
		if (lenID == 0 && !is.null(chipVersion)) {
			## return the entire table based on chipVersion
			mapping <- dbReadTable(conn, chipVersion)
			returnList <- list(chipVersion=chipVersion, species=species, IDType=NA, chipProbeNumber=nrow(mapping), 
					inputProbeNumber=0, matchedProbeNumber=NA)
		} else {
			matchLen <- NULL
			matchField <- NULL
			tableLen <- NULL
			cut0.probeId <- FALSE  # determine whether to cut the 0s in front of the ProbeId. E.g. 0004760445 --> 4760445
			for (i in seq(allTableNames)) {
				tableName.i <- allTableNames[i]
				table.i <- dbReadTable(conn, tableName.i)
				fieldNames.i <- names(table.i)
				## Not include the Symbol column during matching the IDs.
				fieldNames.i <- fieldNames.i[!(fieldNames.i %in% c('Symbol'))]
				len.i <- NULL
				for (fieldNames.ij in fieldNames.i) {
					field.ij <- as.character(table.i[,fieldNames.ij])
					len.ij <- length(which(inputID %in% field.ij))
					if (fieldNames.ij %in% c('ProbeId', 'Array_Address_Id')) {
						len.ij.2 <- length(which(inputID %in% as.character(as.integer(field.ij))))
						if (len.ij.2 > len.ij) {
							cut0.probeId <- TRUE
							len.ij <- len.ij.2
						}
					}
					len.i <- c(len.i, len.ij)
				}
				max.ind.i <- which.max(len.i)
				matchLen <- c(matchLen, len.i[max.ind.i])
				matchField <- c(matchField, fieldNames.i[max.ind.i])
				tableLen <- c(tableLen, nrow(table.i))
			}
			bestMatchLen <- max(matchLen)
			if (bestMatchLen == 0) {
				if (verbose) warning('No matches were found!\n')
				return(c(list(chipVersion=NULL, species=species, IDType=NULL, chipProbeNumber=NULL, 
										matchedProbeNumber=bestMatchLen), idMapping=NULL))
			} else if (bestMatchLen < lenID && verbose) {
				warning('Some input IDs can not be matched!\n')
			}			
			
			if (!returnAllMatches) {
				bestInd <- which(matchLen == bestMatchLen)
				if (length(bestInd) > 1) {
					bestInd <- bestInd[tableLen[bestInd] == min(tableLen[bestInd])]
				}
				matchLen.match <- bestMatchLen
				tableName.match <- allTableNames[bestInd]
				fieldName.match <- matchField[bestInd]
				bestTable <- dbReadTable(conn, tableName.match[1])
				probeNumber <- nrow(bestTable)			
				if (idMapping) {
					## remove duplicated IDs becasue rownames should be unique
					dupInd <- duplicated(bestTable[,fieldName.match[1]])
					bestTable <- bestTable[!dupInd,]
					nuID <- bestTable[,'nuID']
					if (fieldName.match[1] %in% c('ProbeId', 'Array_Address_Id') && cut0.probeId) {
						allID <- as.character(as.integer(bestTable[,fieldName.match[1]]))
					} else {
						allID <- bestTable[,fieldName.match[1]]					
					}
					selInputID <- inputID[inputID %in% allID]
					
					rownames(bestTable) <- allID
					mapping <- as.matrix(bestTable[selInputID, ])
					# if (fieldName.match[1] == 'nuID') {
					# 	bestTable <- bestTable[, names(bestTable) != 'nuID']
					# 	rownames(bestTable) <- nuID
					# 	mapping <- as.matrix(bestTable[selInputID, ])
					# } else {
					# 	names(nuID) <- allID
					# 	mapping <- nuID[selInputID]
					# }
					## match the original input IDs
					mapping <- matchInputID(mapping, inputID.bak)
				}
			} else {
				# if (matchField[which.max(matchLen)] != 'nuID') IDType <- 'nuID'
				matchInd <- which(matchLen > 0)
				tableName.match <- allTableNames[matchInd]
				fieldName.match <- matchField[matchInd]
				matchLen.match <- matchLen[matchInd]
				probeNumber <- NULL
				if (idMapping) mapping <- NULL
				for (i in seq(matchInd)) {
					table.i <- dbReadTable(conn, tableName.match[i])
					if (idMapping) {
						## remove duplicated IDs becasue rownames should be unique
						dupInd.i <- duplicated(table.i[,fieldName.match[i]])
						table.i <- table.i[!dupInd.i,]
						nuID.i <- table.i[,'nuID']
						if (fieldName.match[1] %in% c('ProbeId', 'Array_Address_Id') && cut0.probeId) {
							allID.i <- as.character(as.integer(table.i[,fieldName.match[i]]))
						} else {
							allID.i <- table.i[,fieldName.match[i]]				
						}
						selInputID.i <- inputID[inputID %in% allID.i]
						if (fieldName.match[which.max(matchLen.match)] == 'nuID') {
							table.i <- table.i[, names(table.i) != 'nuID']
							rownames(table.i) <- nuID.i
							mapping.i <- as.matrix(table.i[selInputID.i, ])
							mapping.i <- matchInputID(mapping.i, inputID.bak)
							mapping <- c(mapping, list(mapping.i))
						} else {
							names(nuID.i) <- allID.i
							mapping.i <- nuID[selInputID.i]
							mapping.i <- matchInputID(mapping.i, inputID.bak)
							mapping <- cbind(mapping, mapping.i)
						}
					}
					probeNumber <- c(probeNumber, nrow(table.i))
				}
				rank1 <- rank(matchLen.match)
				rank2 <- rank(-probeNumber)
				ord <- order(rank1 + rank2/length(matchInd), decreasing=TRUE)
				tableName.match <- tableName.match[ord]
				fieldName.match <- fieldName.match[ord]
				probeNumber <- probeNumber[ord]
				matchLen.match <- matchLen.match[ord]
				if (idMapping) {
					if (is(mapping, 'list')) {
						mapping <- mapping[ord]
						names(mapping) <- tableName.match	
					} else if (is(mapping, 'matrix')) {
						rownames(mapping) <- inputID
						mapping <- mapping[,ord]
						colnames(mapping) <- tableName.match		
					}
				}
			}
			returnList <- list(chipVersion=tableName.match, species=species, IDType=fieldName.match, chipProbeNumber=probeNumber, 
					inputProbeNumber=length(inputID), matchedProbeNumber=matchLen.match)
		}
		if (idMapping)  returnList <- c(returnList, list(idMapping=mapping))		
		
		return(returnList)
	}
}

`addNuID2lumi` <-
		function(x.lumi, annotationFile=NULL, sep=NULL, lib.mapping=NULL, annotationColName=c(sequence='Probe_Sequence', target='ILMN_Gene', probe='Probe_Id'), verbose=TRUE) {
	
	history.submitted <- as.character(Sys.time())
	
	## check whether the object is nuID annotated.
	exprs <- exprs(x.lumi)
	id <- rownames(exprs)
	if(all(sapply(id[1:50], is.nuID))) {
		if (verbose) cat('The lumiBatch object is already nuID annotated!\n')
		return(x.lumi)
	}
	if (!is.null(lib.mapping)) {
		if (length(grep('lumi.*\\.db', lib.mapping)) == 0 && length(grep('lumi.*IDMapping', lib.mapping)) == 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!\n'))
			lib.mapping <- NULL
		}
	}
	
	newId <- id
	
	## ---------------------------------------
	## First check whether probe sequence information is available 
	annotation <- pData(featureData(x.lumi))
	names(annotation) <- toupper(names(annotation))
	if (toupper(annotationColName['sequence']) %in% names(annotation)) {
		sequence <- annotation[, toupper(annotationColName['sequence'])]
		cat('Directly converting probe sequence to nuIDs ...\n')
		newId <- sapply(sequence, seq2id)
		names(newId) <- id				
	} else if (!is.null(annotationFile)) {
		## identify the Metadata lines 
		info <- readLines(annotationFile, n=10)    # take the first 10 lines to have a taste
		
		## Use annotationColName[1] as an indicator of Where the metaData stops
		##   intelligently find nMetaDataLines  
		nMetaDataLines <- grep(annotationColName[1], info, ignore.case=TRUE) - 1
		
		if (is.null(sep)) {
			## Find out the separator (sep) by taking the first two line of data, and comparing them.
			##  we assume it is either "," or "\t".
			titleLine <- info[nMetaDataLines + 1]
			sepNum <- regexpr('\t', titleLine)
			if (sepNum >= 2) {
				sep <- '\t'
			} else {
				sepNum <- regexpr(',', titleLine)
				if (sepNum >= 2) {
					sep <- ','
				} else {
					stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!\n')
				}
			}
		}
		
		dataLine1 <- strsplit(info[nMetaDataLines + 2], sep)[[1]]
		quoteCount1 <- gregexpr('"', dataLine1[1])[[1]]
		quoteCount2 <- gregexpr('\'', dataLine1[1])[[1]]
		quote <- ''
		if (sep == ',') quote <- '"'
		if (length(quoteCount1) == 2) {
			quote <- '"'
		} else if (length(quoteCount2) == 2) {
			quote <- '\''
		}
		
		## Read in annotation data
		annotation <- read.table(annotationFile, sep=sep, colClasses="character", header=TRUE, skip=nMetaDataLines,
				blank.lines.skip=TRUE, row.names=NULL, check.names=FALSE, quote=quote, comment.char="", strip.white=TRUE, fill=TRUE)
		
		allColName <- toupper(colnames(annotation))
		colnames(annotation) <- allColName
		## Create unique Id based on 50mer sequence
		if (toupper(annotationColName['sequence']) %in% allColName) {
			nuID <- sapply(annotation[, toupper(annotationColName['sequence'])], seq2id)			
		} else {
			stop('The "sequence" column cannot be found!\nPlease check the "annotationColName" of "sequence"!\n')
		}
		## check the TargetID first
		comm_target <- NULL
		if (toupper(annotationColName['target']) %in% allColName) {
			ann_target <- annotation[, toupper(annotationColName['target'])]
			comm_target <- id[id %in% ann_target]
		}
		if (length(comm_target) == 0) {
			## check the ProbeID if id does not match the TargetID
			if (toupper(annotationColName['probe']) %in% allColName) {
				ann_target <- annotation[, toupper(annotationColName['probe'])]
				comm_target <- id[id %in% ann_target]
				if (length(comm_target) == 0) {
					width <- nchar(ann_target[1])
					id <- formatC(as.numeric(id), width=width, flag='0', format='d')
					comm_target <- id[id %in% ann_target]
					if (length(comm_target) == 0) stop('The annotation file does not match the data!\n')
				}
			}
		} 
		if (length(comm_target) < length(id)) {
			diffId.len <- length(id) - length(comm_target)
			warning(paste('The annotation file does not match the data.',  diffId.len, 'ids cannot be replaced!\n'))
		}
		names(nuID) <- ann_target
		
		newId <- id
		newId[id %in% ann_target] <- nuID[comm_target]
	} else if (!is.null(lib.mapping)) {
		if (!require(lib.mapping, character.only=TRUE)) stop(paste('Annotation library', lib.mapping, 'is not installed!\n'))
		controlId <- c('lysA','pheA','thrB','trpF', 'bla1','bla2','cat1','cat2','cre1','cre2','e1a1',
				'e1a2','gfp1','gfp2','gst1','gst2','gus1','gus2','lux1','lux2','lysA','neo1',
				'neo2','pheA','thrB','trpF')
		if (length(grep('IDMapping', lib.mapping)) == 0) {
			## check the TargetID first
			newId <- mget(id, get(paste(lib.mapping, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
			newId <- unlist(newId)
			if (length(which(!is.na(newId))) == 0) {
				usingTargetID <- FALSE
				## check the ProbeID if id does not match the TargetID
				newId <- mget(id, get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
				if (length(which(!is.na(newId))) == 0) {
					width <- nchar(ls(envir=get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'))[1])
					id <- formatC(as.numeric(id), width=width, flag='0', format='d')
					newId <- mget(id, get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
					if (length(which(!is.na(newId))) == 0) {
						targetID <- pData(featureData(x.lumi))$TargetID
						newId <- mget(targetID, get(paste(lib.mapping, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
						if (length(which(!is.na(newId))) == 0) stop('The library does not match the data!\n')
					}
				} 
			} else {
				usingTargetID <- TRUE
			}
			## Check for the targetIDs cannot be found in the lib.mapping.
			## Some known control genes will not be checked.
			naInd <- which(is.na(newId))
			if (!usingTargetID) {
				TargetID <- featureData(x.lumi)$TargetID
				if (is.null(TargetID)) {
					if (!all(TargetID[naInd] %in% controlId)) {
						if (length(naInd) < 10) {
							warning(paste('Identifiers:', paste(TargetID[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
						} else {
							warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
						}
					}
				}
			} else if (!all(id[naInd] %in% controlId)) {
				if (length(naInd) < 10) {
					warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
				} else {
					warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
				}
			}
			newId[naInd] <- id[naInd]
		} else {
			chipInfo <- getChipInfo(featureNames(x.lumi), lib.mapping=lib.mapping, idMapping=TRUE, verbose=FALSE)
			newId <- chipInfo$idMapping[,'nuID']
			naInd <- which(is.na(newId))
			if (length(naInd) > 0) {
				newId[naInd] <- id[naInd]
				if (!all(id[naInd] %in% controlId))	{
					if (length(naInd) > 500) {
						warning(paste('More than 500 identifiers cannot be found in the library:', lib.mapping, 
										'!\n The provided library might be wrong or outdated!\n', sep=''))
					} else if (length(naInd) < 10) {
						warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
					} else {
						warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
					}
				}
			}
			selInfoName <- names(chipInfo)
			selInfoName <- selInfoName[selInfoName != 'idMapping']
			chipInfo.print <- paste(selInfoName, unlist(chipInfo[selInfoName]), sep=': ')
			notes(x.lumi) <- c(notes(x.lumi), list('Chip Information'=chipInfo.print))
			
			annotation(x.lumi) <- switch(chipInfo$species,
					'Human'='lumiHumanAll.db',
					'Mouse'='lumiMouseAll.db',
					'Rat'='lumiRatAll.db')
		}
	} else {
		cat('Please provide the annotation file or lumi annotation library!\n')
	}
	if (all(newId == id)) {
		conversion <- FALSE
	} else {
		conversion <- TRUE
	}
	
	if (any(duplicated(newId)))  {
		if (verbose) cat('Duplicated IDs found and were merged!\n')
		dupId <- unique(newId[duplicated(newId)])
		## determine whether the detection p-value close to 0 or 1 is significant
		if (!is.null(detection(x.lumi))) {
			detect.low <- exprs[which.max(detection(x.lumi)[,1]), 1]
			detect.high <- exprs[which.min(detection(x.lumi)[,1]), 1]
		}
		
		rmIndex <- NULL
		for (dupId.i in dupId) {
			dupIndex <- which(newId == dupId.i)
			ave.exp <- colMeans(exprs(x.lumi)[dupIndex,, drop=FALSE])
			exprs(x.lumi)[dupIndex[1],] <- ave.exp
			if (is(x.lumi, 'LumiBatch') && !is.null(beadNum(x.lumi)) && !is.null(detection(x.lumi))) {
				totalBeadNum <- colSums(beadNum(x.lumi)[dupIndex, , drop=FALSE])
				if (detect.low < detect.high) {
					maxDetection <- apply(detection(x.lumi), 2, min)
				} else {
					maxDetection <- apply(detection(x.lumi), 2, max)
				}
				
				temp <- colSums(se.exprs(x.lumi)[dupIndex, , drop=FALSE]^2 * (beadNum(x.lumi)[dupIndex,, drop=FALSE] - 1))
				temp <- temp / (totalBeadNum - length(dupIndex))
				se.exprs(x.lumi)[dupIndex[1],] <- sqrt(temp * (colSums(1/beadNum(x.lumi)[dupIndex,, drop=FALSE])))
				detection(x.lumi)[dupIndex[1],] <- maxDetection
				beadNum(x.lumi)[dupIndex[1],] <- totalBeadNum
			}
			rmIndex <- c(rmIndex, dupIndex[-1])
		}
		## remove duplicated
		x.lumi <- x.lumi[-rmIndex, ]
		newId <- newId[-rmIndex]
	}
	
	## update the feature names (probe ids)
	if (conversion) {
		names(newId) <- NULL
		featureNames(x.lumi) <- newId
		## update the feautre data
		featureData <- featureData(x.lumi)
		rownames(pData(featureData)) <- newId
		if (!is.null(pData(featureData)$'PROBE_SEQUENCE')) pData(featureData)$'PROBE_SEQUENCE' <- NULL
		featureData(x.lumi) <- featureData
		
		## Add history tracking
		if (is(x.lumi, 'LumiBatch')) {
			history.finished <- as.character(Sys.time())
			history.command <- capture.output(print(match.call(addNuID2lumi)))
			if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
							finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	}
	
	return(x.lumi)
}

`addNuId2lumi` <-
		function(x.lumi, annotationFile=NULL, sep=NULL, lib.mapping=NULL, annotationColName=c(sequence='Probe_Sequence', target='Target', probe='Probe_Id'), verbose=TRUE) {
	cat('Function addNuId2lumi is deprecated!\n Please use addNuID2lumi instead!\n')
	addNuID2lumi(x.lumi, annotationFile, sep, lib.mapping, annotationColName, verbose)
}

x.lumi<-addNuID2lumi(x.lumi, lib.mapping=mapping)

# Quality control (not run)
# QC
# q.lumi <- lumiQ(x.lumi)  

# Background correction
if(background.correction=="bgAdjust.affy") {
   x.lumi<-lumiB(x.lumi, method="bgAdjust.affy")
}

# Transformation
if(transformation=="none") {
   t.lumi<-x.lumi
}
if(transformation=="vst") {
   t.lumi<-lumiT(x.lumi, method=c("vst"))
}
if(transformation=="log2") {
   t.lumi<-lumiT(x.lumi, method=c("log2"))
}

# Normalization
if(normalize.chips=="none") {
   n.lumi<-t.lumi
}
if(normalize.chips=="rsn") {
   n.lumi<-lumiN(t.lumi, method=c("rsn"))
}
if(normalize.chips=="loess") {
   n.lumi<-lumiN(t.lumi, method=c("loess"))
}
if(normalize.chips=="quantile") {
   n.lumi<-lumiN(t.lumi, method=c("quantile"))
}
if(normalize.chips=="vsn") {
   n.lumi<-lumiN(t.lumi, method=c("vsn"))
}

# Convert sample names to Chipster style
dat2<-exprs(n.lumi)
sample.names<-colnames(dat2)
sample.names<-paste("chip.", sample.names, sep="")
names(dat2)<-sample.names
colnames(dat2)<-sample.names

# Write out a phenodata
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
write.table(data.frame(sample=sample.names, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

if(chiptype!="Illumina") {
   # Including gene names to data
   library(chiptype, character.only=T)
   # symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "SYMBOL", sep="")))))[x.lumi@featureData@data$TargetID,])
   # genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "GENENAME", sep="")))))[x.lumi@featureData@data$TargetID,])
   # symbol<-gsub("#", "", symbol)
   # genename<-gsub("#", "", genename)
   # Write out expression data
   # write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	symbols <- unlist (lookUp(rownames(dat2), chiptype, what="SYMBOL"))
	genenames <- unlist (lookUp(rownames(dat2), chiptype, what="GENENAME"))
	symbols <- gsub("#", "", symbols)
	genenames <- gsub("#", "", genenames)
	symbols <- gsub("'", "", symbols)
	genenames <- gsub("'", "", genenames)
	write.table(data.frame(symbol=symbols, description=genenames, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} else {
   # Write out expression data
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
