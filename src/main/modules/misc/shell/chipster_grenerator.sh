#!/bin/bash

tool=$(basename $1 .acd)

java -cp chipster-all-2.9.1.jar fi.csc.microarray.analyser.emboss.ACDToSADL $tool.acd $tool.R
awk '{print "\# "$0}' $tool.R > sadl

rm -f  $tool.R

echo 'emboss.parameters <- paste("-auto")' > emboss.param 
sed s/"OPTIONAL"/""/g sadl | sed s/:/""/g | awk '{print "emboss.parameters \<\- paste\(emboss.parameters\, \"\-"$3"\", "$3 ")"}' >> emboss.param 

echo '# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)' >> sadl



cat <<EOF  > alku

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
	stop(paste('Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "$tool")
EOF


cat <<EOF > loppu

command.full <- paste(emboss.binary, emboss.parameters, ' >> $tool.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> $tool.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f $tool.log")
}

EOF

cat sadl alku emboss.param loppu >  $tool.R

rm sadl alku emboss.param loppu
