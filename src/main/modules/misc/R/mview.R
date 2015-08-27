# TOOL mview.R: "Display sequence alignment or BLAST report as HTML page" (Display sequence alignment or BLAST report as a coloured HTML page)
# INPUT sequences: sequences TYPE GENERIC
# OUTPUT OPTIONAL alignment.html
# OUTPUT OPTIONAL mview.log
# PARAMETER moltype: "Molecule type" TYPE [aa: Protein, dna: DNA, rna: RNA, na: NA] DEFAULT aa (Molecule type. Affects coloring.)
# PARAMETER inputtype: "Alignment format" TYPE [clustal: Clustal, fasta: Fasta, blast: "BLAST report"] DEFAULT clustal (The format of input alignment.)
# PARAMETER OPTIONAL width: "Row length" TYPE STRING DEFAULT "flat" (The define the row length in the alignment. Flat: infine row length.) 
# PARAMETER OPTIONAL ruler: "Show ruler" TYPE [off: No, on: Yes] DEFAULT off (Show ruler.)
# PARAMETER OPTIONAL conservation: "Show conservation line" TYPE [off: No, on: Yes] DEFAULT off (Show clustal conservation line.)
# PARAMETER OPTIONAL consensus: "Show consensus" TYPE [off: No, on: Yes] DEFAULT off (Show consensus sequence.)
# PARAMETER OPTIONAL pcid: "Identity precent calculation method." TYPE [aligned: Aligned, hit: Hit, reference: Reference] DEFAULT aligned (Select the reference method for computing percent identities.)
# PARAMETER OPTIONAL reference: "Select reference sequence" TYPE STRING DEFAULT "query" (Use row number or row identifier as identity precent calculation reference.)
# PARAMETER OPTIONAL coloring: "Coloring style" TYPE [any: All, consensus: "Consensus only", groups: Groups, identity: "Identical columns only", none: None] DEFAULT any (Define which parts of the alignment will be colored)
# PARAMETER OPTIONAL colormap: "Color map to use" TYPE [CCLUSTAL: "highlight equivalence class", CHARGE: "highlight charged amino acids", CLUSTAL: "highlight amino acid physicochemical properties", CLUSTAL_NUC: "CLUSTAL-derived colours for nucleotides", CYS: "Highlight cysteines", D1: "Highlight nucleotide types", D2: "DNA: highlight match versus mismatch" ] DEFAULT CLUSTAL (Colormap to use)
# PARAMETER OPTIONAL css: "Color background" TYPE [on: Yes, off: no] DEFAULT off (Color the backgroud of sequeces in stead of letters )
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 20.07. 2015
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
mview.path <- file.path(chipster.tools.path, "mview" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequences")

#check sequece file type
if (inputtype == "clustal" || inputtype == "fasta" ){
  inputfile.to.check <- ("sequences")
  sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
  sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
  str.filetype <- system(sfcheck.command, intern = TRUE )

  if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	  stop("CHIPSTER-NOTE: The input file does not match the selected input format.")
  }

  #count the query sequeces
  seqcount.exe <- file.path(emboss.path, "seqcount sequences -filter")
  str.queryseq <- system(seqcount.exe, intern = TRUE )
  num.queryseq <- as.integer(str.queryseq)

  #round(num.queryseq)

  if (num.queryseq > 1000){
	stop(paste('CHIPSTER-NOTE: Too many sequences in the alignmnet. Maximun is 1000 but your file contains ', num.queryseq ))
  }
}

if (inputtype == "blast" ){
	nohtml.binary <- file.path(emboss.path, "nohtml")
	system("mv sequences sequences_blast ")
	nohtml.command <- paste(nohtml.binary, "sequences_blast", "sequences")
	system(nohtml.command)
}

  
mview.binary <- file.path(mview.path, "mview")
mview.parameters <- paste("-html full -bold")
mview.parameters <- paste(mview.parameters, "-in", inputtype)
mview.parameters <- paste(mview.parameters, "-ruler", ruler)
mview.parameters <- paste(mview.parameters, "-pcid", pcid)
mview.parameters <- paste(mview.parameters, "-reference", reference)
mview.parameters <- paste(mview.parameters, "-conservation", conservation)
mview.parameters <- paste(mview.parameters, "-consensus", consensus)
mview.parameters <- paste(mview.parameters, "-coloring", coloring)
mview.parameters <- paste(mview.parameters, "-colormap", colormap)
mview.parameters <- paste(mview.parameters, "-width", width)
mview.parameters <- paste(mview.parameters, "-css", css)	



command.full <- paste(mview.binary, mview.parameters, ' < sequences > alignment.html 2>> mview.log' )
echo.command <- paste('echo "',command.full, ' "> mview.log' )
system(echo.command)

system(command.full)
system("ls -l >> mview.log")

if ( save_log == "no") {
	system ("rm -f mview.log")
}