# TOOL infoalign.R: "Alignment statistics" (Display basic information about a multiple sequence alignment)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT alignment_summary.html
# OUTPUT OPTIONAL infoalign.log
# PARAMETER OPTIONAL refseq: "The number or the name of the reference sequence" TYPE STRING DEFAULT 0 (If you give the number in the alignment or the name of a sequence, it will be taken to be the reference sequence. The reference sequence is the one against which all the other sequences are compared. If this is set to 0 then the consensus sequence will be used as the reference sequence. By default the consensus sequence is used as the reference sequence.)
# PARAMETER OPTIONAL plurality: "Plurality check % for consensus" TYPE DECIMAL FROM 0.0 TO 100.0 DEFAULT 50.0 (Set a cut-off for the % of positive scoring matches below which there is no consensus. The default plurality is taken as 50% of the total weight of all the sequences in the alignment.)
# PARAMETER OPTIONAL identity: "Required % of identities at a position for consensus" TYPE DECIMAL FROM 0.0 TO 100.0 DEFAULT 0.0 (Provides the facility of setting the required number of identities at a position for it to give a consensus. Therefore, if this is set to 100% only columns of identities contribute to the consensus.)
# PARAMETER OPTIONAL heading: "Display column headings" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display column headings)
# PARAMETER OPTIONAL name: "Display 'name' column" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display 'name' column)
# PARAMETER OPTIONAL seqlength: "Display 'seqlength' column" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display 'seqlength' column)
# PARAMETER OPTIONAL alignlength: "Display 'alignlength' column" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display 'alignlength' column)
# PARAMETER OPTIONAL gaps: "Display number of gaps" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display number of gaps)
# PARAMETER OPTIONAL gapcount: "Display number of gap positions" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display number of gap positions)
# PARAMETER OPTIONAL idcount: "Display number of identical positions" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display number of identical positions)
# PARAMETER OPTIONAL simcount: "Display number of similar positions" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display number of similar positions)
# PARAMETER OPTIONAL diffcount: "Display number of different positions" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display number of different positions)
# PARAMETER OPTIONAL change: "Display % number of changed positions" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display % number of changed positions)
# PARAMETER OPTIONAL weight: "Display 'weight' column" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display 'weight' column)
# PARAMETER OPTIONAL description: "Display 'description' column" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display 'description' column)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("sequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "infoalign")
emboss.parameters <- paste("-auto -html Y -only")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outfile alignment_summary.html")
emboss.parameters <- paste(emboss.parameters, "-refseq", refseq)
emboss.parameters <- paste(emboss.parameters, "-plurality", plurality)
emboss.parameters <- paste(emboss.parameters, "-identity", identity)
emboss.parameters <- paste(emboss.parameters, "-heading", heading)
emboss.parameters <- paste(emboss.parameters, "-name", name)
emboss.parameters <- paste(emboss.parameters, "-seqlength", seqlength)
emboss.parameters <- paste(emboss.parameters, "-alignlength", alignlength)
emboss.parameters <- paste(emboss.parameters, "-gaps", gaps)
emboss.parameters <- paste(emboss.parameters, "-gapcount", gapcount)
emboss.parameters <- paste(emboss.parameters, "-idcount", idcount)
emboss.parameters <- paste(emboss.parameters, "-simcount", simcount)
emboss.parameters <- paste(emboss.parameters, "-diffcount", diffcount)
emboss.parameters <- paste(emboss.parameters, "-change", change)
emboss.parameters <- paste(emboss.parameters, "-weight", weight)
emboss.parameters <- paste(emboss.parameters, "-description", description)

command.full <- paste(emboss.binary, emboss.parameters, ' >> infoalign.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> infoalign.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f infoalign.log")
}
