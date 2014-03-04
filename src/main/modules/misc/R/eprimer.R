# TOOL eprimer.R: "Design PCR primers and hybridization oligos with Primer3" (Picks PCR primers and hybridization oligos)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL primers.txt
# OUTPUT OPTIONAL eprimer3.log 
# PARAMETER OPTIONAL primer: "Pick PCR primer(s\)" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Tell EPrimer3 to pick primer(s\))
# PARAMETER OPTIONAL task: "Select task" TYPE [1: "Pick PCR primers", 2: "Pick forward primer only", 3: "Pick reverse primer only", 4: "No primers needed"] FROM 1 TO 1 DEFAULT 1 (Tell EPrimer3 what task to perform. Legal values are 1: 'Pick PCR primers', 2: 'Pick forward primer only', 3: 'Pick reverse primer only', 4: 'No primers needed'.)
# PARAMETER OPTIONAL hybridprobe: "Pick hybridization probe" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (An 'internal oligo' is intended to be used as a hybridization probe (hyb probe\) to detect the PCR product after amplification.)
# PARAMETER OPTIONAL numreturn: "Number of results to return" TYPE INTEGER FROM 0 DEFAULT 5 (The maximum number of primer pairs to return. Primer pairs returned are sorted by their 'quality', in other words by the value of the objective function (where a lower number indicates a better primer pair\). Caution: setting this parameter to a large value will increase running time.)
# PARAMETER OPTIONAL includedregion: "Included region(s\)" TYPE STRING (A sub-region of the given sequence in which to pick primers. For example, often the first dozen or so bases of a sequence are vector, and should be excluded from consideration. The value for this parameter has the form \ (start\),(end\) \ where (start\) is the index of the first base to consider, and (end\) is the last in the primer-picking region.)
# PARAMETER OPTIONAL targetregion: "Target region(s\)" TYPE STRING (If one or more Targets is specified then a legal primer pair must flank at least one of them. A Target might be a simple sequence repeat site (for example a CA repeat\) or a single-base-pair polymorphism. The value should be a space-separated list of \ (start\),(end\) \ pairs where (start\) is the index of the first base of a Target, and (end\) is the last \ E.g. 50,51 requires primers to surround the 2 bases at positions 50 and 51.)
# PARAMETER OPTIONAL excludedregion: "Excluded region(s\)" TYPE STRING (Primer oligos may not overlap any region specified in this tag. The associated value must be a space-separated list of \ (start\),(end\) \ pairs where (start\) is the index of the first base of the excluded region, and and (end\) is the last. This tag is useful for tasks such as excluding regions of low sequence quality or for excluding regions containing repetitive elements such as ALUs or LINEs. \ E.g. 401,407 68,70 forbids selection of primers in the 7 bases starting at 401 and the 3 bases at 68.)
# PARAMETER OPTIONAL forwardinput: "Forward input primer sequence to check" TYPE STRING (The sequence of a forward primer to check and around which to design reverse primers and optional internal oligos. Must be a substring of SEQUENCE.)
# PARAMETER OPTIONAL reverseinput: "Reverse input primer sequence to check" TYPE STRING (The sequence of a reverse primer to check and around which to design forward primers and optional internal oligos. Must be a substring of the reverse strand of SEQUENCE.)
# PARAMETER OPTIONAL gcclamp: "GC clamp" TYPE INTEGER FROM 0 DEFAULT 0 (Require the specified number of consecutive Gs and Cs at the 3' end of both the forward and reverse primer. (This parameter has no effect on the internal oligo if one is requested.\))
# PARAMETER OPTIONAL osize: "Primer optimum size" TYPE INTEGER FROM 0 DEFAULT 20 (Optimum length (in bases\) of a primer oligo. EPrimer3 will attempt to pick primers close to this length.)
# PARAMETER OPTIONAL minsize: "Primer minimum size" TYPE INTEGER FROM 1 DEFAULT 18 (Minimum acceptable length of a primer. Must be greater than 0 and less than or equal to MAX-SIZE.)
# PARAMETER OPTIONAL maxsize: "Primer maximum size" TYPE INTEGER FROM 18 TO 99 DEFAULT 27 (Maximum acceptable length (in bases\) of a primer. Currently this parameter should not be larger than 35. This limit is governed by the maximum oligo size for which EPrimer3's melting-temperature is valid.)
# PARAMETER OPTIONAL opttm: "Primer optimum Tm" TYPE DECIMAL DEFAULT 60.0 (Optimum melting temperature(Celsius\) for a primer oligo. EPrimer3 will try to pick primers with melting temperatures are close to this temperature. The oligo melting temperature formula in EPrimer3 is that given in Rychlik, Spencer and Rhoads, Nucleic Acids Research, vol 18, num 21, pp 6409-6412 and Breslauer, Frank, Bloecker and Marky, Proc. Natl. Acad. Sci. USA, vol 83, pp 3746-3750. Please refer to the former paper for background discussion.)
# PARAMETER OPTIONAL mintm: "Primer minimum Tm" TYPE DECIMAL DEFAULT 57.0 (Minimum acceptable melting temperature(Celsius\) for a primer oligo.)
# PARAMETER OPTIONAL maxtm: "Primer maximum Tm" TYPE DECIMAL FROM 57.0 DEFAULT 63.0 (Maximum acceptable melting temperature(Celsius\) for a primer oligo.)
# PARAMETER OPTIONAL maxdifftm: "Maximum difference in Tm of primers" TYPE DECIMAL DEFAULT 100.0 (Maximum acceptable (unsigned\) difference between the melting temperatures of the forward and reverse primers.)
# PARAMETER OPTIONAL ogcpercent: "Primer optimum GC percent" TYPE DECIMAL DEFAULT 50.0 (Primer optimum GC percent.)
# PARAMETER OPTIONAL mingc: "Primer minimum GC percent" TYPE DECIMAL DEFAULT 20.0 (Minimum allowable percentage of Gs and Cs in any primer.)
# PARAMETER OPTIONAL maxgc: "Primer maximum GC percent" TYPE DECIMAL FROM 20.0 DEFAULT 80.0 (Maximum allowable percentage of Gs and Cs in any primer generated by Primer.)
# PARAMETER OPTIONAL saltconc: "Salt concentration (mM\)" TYPE DECIMAL DEFAULT 50.0 (The millimolar concentration of salt (usually KCl\) in the PCR. EPrimer3 uses this argument to calculate oligo melting temperatures.)
# PARAMETER OPTIONAL dnaconc: "DNA concentration (nM\)" TYPE DECIMAL DEFAULT 50.0 (The nanomolar concentration of annealing oligos in the PCR. EPrimer3 uses this argument to calculate oligo melting temperatures. The default (50nM\) works well with the standard protocol used at the Whitehead/MIT Center for Genome Research--0.5 microliters of 20 micromolar concentration for each primer oligo in a 20 microliter reaction with 10 nanograms template, 0.025 units/microliter Taq polymerase in 0.1 mM each dNTP, 1.5mM MgCl2, 50mM KCl, 10mM Tris-HCL (pH 9.3\) using 35 cycles with an annealing temperature of 56 degrees Celsius. This parameter corresponds to 'c' in Rychlik, Spencer and Rhoads' equation (ii\) (Nucleic Acids Research, vol 18, num 21\) where a suitable value (for a lower initial concentration of template\) is 'empirically determined'. The value of this parameter is less than the actual concentration of oligos in the reaction because it is the concentration of annealing oligos, which in turn depends on the amount of template (including PCR product\) in a given cycle. This concentration increases a great deal during a PCR; fortunately PCR seems quite robust for a variety of oligo melting temperatures. \ See ADVICE FOR PICKING PRIMERS.)
# PARAMETER OPTIONAL maxpolyx: "Maximum polynucleotide repeat" TYPE INTEGER FROM 0 DEFAULT 5 (The maximum allowable length of a mononucleotide repeat in a primer, for example AAAAAA.)
# PARAMETER OPTIONAL psizeopt: "Product optimum size" TYPE INTEGER FROM 0 DEFAULT 200 (The optimum size for the PCR product. 0 indicates that there is no optimum product size.)
# PARAMETER OPTIONAL prange: "Product size range" TYPE STRING DEFAULT 100-300 (The associated values specify the lengths of the product that the user wants the primers to create, and is a space separated list of elements of the form \ (x\)-(y\) \ where an (x\)-(y\) pair is a legal range of lengths for the product. For example, if one wants PCR products to be between 100 to 150 bases (inclusive\) then one would set this parameter to 100-150. If one desires PCR products in either the range from 100 to 150 bases or in the range from 200 to 250 bases then one would set this parameter to 100-150 200-250. \ EPrimer3 favors ranges to the left side of the parameter string. EPrimer3 will return legal primers pairs in the first range regardless the value of the objective function for these pairs. Only if there are an insufficient number of primers in the first range will EPrimer3 return primers in a subsequent range.)
# PARAMETER OPTIONAL ptmopt: "Product optimum Tm" TYPE DECIMAL DEFAULT 0.0 (The optimum melting temperature for the PCR product. 0 indicates that there is no optimum temperature.)
# PARAMETER OPTIONAL ptmmin: "Product minimum Tm" TYPE DECIMAL DEFAULT -1000000.0 (The minimum allowed melting temperature of the amplicon. Please see the documentation on the maximum melting temperature of the product for details.)
# PARAMETER OPTIONAL ptmmax: "Product maximum Tm" TYPE DECIMAL FROM -1000000.0 DEFAULT 1000000.0 (The maximum allowed melting temperature of the amplicon. Product Tm is calculated using the formula from Bolton and McCarthy, PNAS 84:1390 (1962\) as presented in Sambrook, Fritsch and Maniatis, Molecular Cloning, p 11.46 (1989, CSHL Press\). \ Tm = 81.5 + 16.6(log10([Na+]\)\) + .41*(%GC\) - 600/length \ Where [Na+} is the molar sodium concentration, (%GC\) is the percent of Gs and Cs in the sequence, and length is the length of the sequence. \ A similar formula is used by the prime primer selection program in GCG, which instead uses 675.0/length in the last term (after F. Baldino, Jr, M.-F. Chesselet, and M.E. Lewis, Methods in Enzymology 168:766 (1989\) eqn (1\) on page 766 without the mismatch and formamide terms\). The formulas here and in Baldino et al. assume Na+ rather than K+. According to J.G. Wetmur, Critical Reviews in BioChem. and Mol. Bio. 26:227 (1991\) 50 mM K+ should be equivalent in these formulae to .2 M Na+. EPrimer3 uses the same salt concentration value for calculating both the primer melting temperature and the oligo melting temperature. If you are planning to use the PCR product for hybridization later this behavior will not give you the Tm under hybridization conditions.)
# PARAMETER OPTIONAL oexcludedregion: "Internal oligo excluded region" TYPE STRING (Middle oligos may not overlap any region specified by this tag. The associated value must be a space-separated list of \ (start\),(end\) \ pairs, where (start\) is the index of the first base of an excluded region, and (end\) is the last. Often one would make Target regions excluded regions for internal oligos.)
# PARAMETER OPTIONAL oligoinput: "Internal oligo input sequence (if any\)" TYPE STRING (The sequence of an internal oligo to check and around which to design forward and reverse primers. Must be a substring of SEQUENCE.)
# PARAMETER OPTIONAL osizeopt: "Internal oligo optimum size" TYPE INTEGER FROM 0 DEFAULT 20 (Optimum length (in bases\) of an internal oligo. EPrimer3 will attempt to pick primers close to this length.)
# PARAMETER OPTIONAL ominsize: "Internal oligo minimum size" TYPE INTEGER FROM 0 DEFAULT 18 (Minimum acceptable length of an internal oligo. Must be greater than 0 and less than or equal to INTERNAL-OLIGO-MAX-SIZE.)
# PARAMETER OPTIONAL omaxsize: "Internal oligo maximum size" TYPE INTEGER FROM 18 TO 35 DEFAULT 27 (Maximum acceptable length (in bases\) of an internal oligo. Currently this parameter cannot be larger than 35. This limit is governed by maximum oligo size for which EPrimer3's melting-temperature is valid.)
# PARAMETER OPTIONAL otmopt: "Internal oligo optimum Tm" TYPE DECIMAL DEFAULT 60.0 (Optimum melting temperature (Celsius\) for an internal oligo. EPrimer3 will try to pick oligos with melting temperatures that are close to this temperature. The oligo melting temperature formula in EPrimer3 is that given in Rychlik, Spencer and Rhoads, Nucleic Acids Research, vol 18, num 21, pp 6409-6412 and Breslauer, Frank, Bloecker and Marky, Proc. Natl. Acad. Sci. USA, vol 83, pp 3746-3750. Please refer to the former paper for background discussion.)
# PARAMETER OPTIONAL otmmin: "Internal oligo minimum Tm" TYPE DECIMAL DEFAULT 57.0 (Minimum acceptable melting temperature(Celsius\) for an internal oligo.)
# PARAMETER OPTIONAL otmmax: "Internal oligo maximum Tm" TYPE DECIMAL FROM 57.0 DEFAULT 63.0 (Maximum acceptable melting temperature (Celsius\) for an internal oligo.)
# PARAMETER OPTIONAL ogcopt: "Internal oligo optimum GC percent" TYPE DECIMAL DEFAULT 50.0 (Internal oligo optimum GC percent.)
# PARAMETER OPTIONAL ogcmin: "Internal oligo minimum GC" TYPE DECIMAL DEFAULT 20.0 (Minimum allowable percentage of Gs and Cs in an internal oligo.)
# PARAMETER OPTIONAL ogcmax: "Internal oligo maximum GC" TYPE DECIMAL FROM 20.0 DEFAULT 80.0 (Maximum allowable percentage of Gs and Cs in any internal oligo generated by Primer.)
# PARAMETER OPTIONAL osaltconc: "Internal oligo salt concentration (mM\)" TYPE DECIMAL DEFAULT 50.0 (The millimolar concentration of salt (usually KCl\) in the hybridization. EPrimer3 uses this argument to calculate internal oligo melting temperatures.)
# PARAMETER OPTIONAL odnaconc: "Internal oligo DNA concentration (nM\)" TYPE DECIMAL DEFAULT 50.0 (The nanomolar concentration of annealing internal oligo in the hybridization.)
# PARAMETER OPTIONAL oanyself: "Internal oligo maximum self complementarity" TYPE DECIMAL TO 9999.99 DEFAULT 12.00 (The maximum allowable local alignment score when testing an internal oligo for (local\) self-complementarity. Local self-complementarity is taken to predict the tendency of oligos to anneal to themselves The scoring system gives 1.00 for complementary bases, -0.25 for a match of any base (or N\) with an N, -1.00 for a mismatch, and -2.00 for a gap. Only single-base-pair gaps are allowed. For example, the alignment \ 5' ATCGNA 3' \ || | | \ 3' TA-CGT 5' \ is allowed (and yields a score of 1.75\), but the alignment \ 5' ATCCGNA 3' \ || | | \ 3' TA--CGT 5' \ is not considered. Scores are non-negative, and a score of 0.00 indicates that there is no reasonable local alignment between two oligos.)
# PARAMETER OPTIONAL oendself: "Internal oligo maximum 3' self complementarity" TYPE DECIMAL FROM 12.00 TO 9999.99 DEFAULT 12.00 (The maximum allowable 3'-anchored global alignment score when testing a single oligo for self-complementarity. \ The scoring system is as for the Maximum Complementarity argument. In the examples above the scores are 7.00 and 6.00 respectively. Scores are non-negative, and a score of 0.00 indicates that there is no reasonable 3'-anchored global alignment between two oligos. In order to estimate 3'-anchored global alignments for candidate oligos, Primer assumes that the sequence from which to choose oligos is presented 5' to 3'. \ INTERNAL-OLIGO-SELF-END is meaningless when applied to internal oligos used for hybridization-based detection, since primer-dimer will not occur. We recommend that INTERNAL-OLIGO-SELF-END be set at least as high as INTERNAL-OLIGO-SELF-ANY.)
# PARAMETER OPTIONAL opolyxmax: "Internal oligo maximum polynucleotide repeat" TYPE INTEGER FROM 0 DEFAULT 5 (The maximum allowable length of an internal oligo mononucleotide repeat, for example AAAAAA.)
# PARAMETER OPTIONAL omishybmax: "Internal oligo maximum mishybridization" TYPE DECIMAL TO 9999.99 DEFAULT 12.0 (Similar to MAX-MISPRIMING except that this parameter applies to the similarity of candidate internal oligos to the library specified in INTERNAL-OLIGO-MISHYB-LIBRARY.)
# PARAMETER OPTIONAL explainflag: "Explain flag" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this flag is true, produce LEFT-EXPLAIN, RIGHT-EXPLAIN, and INTERNAL-OLIGO-EXPLAIN output tags, which are intended to provide information on the number of oligos and primer pairs that EPrimer3 examined, and statistics on the number discarded for various reasons.)
# PARAMETER OPTIONAL fileflag: "Create results files for each sequence" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If the associated value is true, then EPrimer3 creates two output files for each input SEQUENCE. File (sequence-id\).for lists all acceptable forward primers for (sequence-id\), and (sequence-id\).rev lists all acceptable reverse primers for (sequence-id\), where (sequence-id\) is the value of the SEQUENCE-ID tag (which must be supplied\). In addition, if the input tag TASK is 1 or 4, EPrimer3 produces a file (sequence-id\).int, which lists all acceptable internal oligos.)
# PARAMETER OPTIONAL firstbaseindex: "First base index" TYPE INTEGER DEFAULT 1 (This parameter is the index of the first base in the input sequence. For input and output using 1-based indexing (such as that used in GenBank and to which many users are accustomed\) set this parameter to 1. For input and output using 0-based indexing set this parameter to 0. (This parameter also affects the indexes in the contents of the files produced when the primer file flag is set.\))
# PARAMETER OPTIONAL pickanyway: "Pick anyway" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If true pick a primer pair even if LEFT-INPUT, RIGHT-INPUT, or INTERNAL-OLIGO-INPUT violates specific constraints.)
# PARAMETER OPTIONAL maxmispriming: "Primer maximum mispriming" TYPE DECIMAL TO 9999.99 DEFAULT 12.00 (The maximum allowed weighted similarity with any sequence in MISPRIMING-LIBRARY.)
# PARAMETER OPTIONAL pairmaxmispriming: "Primer pair maximum mispriming" TYPE DECIMAL TO 9999.99 DEFAULT 24.00 (The maximum allowed sum of weighted similarities of a primer pair (one similarity for each primer\) with any single sequence in MISPRIMING-LIBRARY.)
# PARAMETER OPTIONAL numnsaccepted: "Maximum Ns accepted in a primer" TYPE INTEGER FROM 0 DEFAULT 0 (Maximum number of unknown bases (N\) allowable in any primer.)
# PARAMETER OPTIONAL selfany: "Maximum self complementarity" TYPE DECIMAL FROM 0.00 TO 9999.99 DEFAULT 8.00 (The maximum allowable local alignment score when testing a single primer for (local\) self-complementarity and the maximum allowable local alignment score when testing for complementarity between forward and reverse primers. Local self-complementarity is taken to predict the tendency of primers to anneal to each other without necessarily causing self-priming in the PCR. The scoring system gives 1.00 for complementary bases, -0.25 for a match of any base (or N\) with an N, -1.00 for a mismatch, and -2.00 for a gap. Only single-base-pair gaps are allowed. For example, the alignment \ 5' ATCGNA 3' \ ...|| | | \ 3' TA-CGT 5' \ is allowed (and yields a score of 1.75\), but the alignment \ 5' ATCCGNA 3' \ ...|| | | \ 3' TA--CGT 5' \ is not considered. Scores are non-negative, and a score of 0.00 indicates that there is no reasonable local alignment between two oligos.)
# PARAMETER OPTIONAL selfend: "Maximum 3' self complementarity" TYPE DECIMAL FROM 0.00 TO 8.00 DEFAULT 3.00 (The maximum allowable 3'-anchored global alignment score when testing a single primer for self-complementarity, and the maximum allowable 3'-anchored global alignment score when testing for complementarity between forward and reverse primers. The 3'-anchored global alignment score is taken to predict the likelihood of PCR-priming primer-dimers, for example \ 5' ATGCCCTAGCTTCCGGATG 3' \ .............||| ||||| \ ..........3' AAGTCCTACATTTAGCCTAGT 5' \ or \ 5' AGGCTATGGGCCTCGCGA 3' \ ...............|||||| \ ............3' AGCGCTCCGGGTATCGGA 5' \ The scoring system is as for the Maximum Complementarity argument. In the examples above the scores are 7.00 and 6.00 respectively. Scores are non-negative, and a score of 0.00 indicates that there is no reasonable 3'-anchored global alignment between two oligos. In order to estimate 3'-anchored global alignments for candidate primers and primer pairs, Primer assumes that the sequence from which to choose primers is presented 5' to 3'. It is nonsensical to provide a larger value for this parameter than for the Maximum (local\) Complementarity parameter because the score of a local alignment will always be at least as great as the score of a global alignment.)
# PARAMETER OPTIONAL maxendstability: "Maximum 3' end stability" TYPE DECIMAL TO 999.9999 DEFAULT 9.0 (The maximum stability for the five 3' bases of a forward or reverse primer. Bigger numbers mean more stable 3' ends. The value is the maximum delta G for duplex disruption for the five 3' bases as calculated using the nearest neighbor parameters published in Breslauer, Frank, Bloecker and Marky, Proc. Natl. Acad. Sci. USA, vol 83, pp 3746-3750. EPrimer3 uses a completely permissive default value for backward compatibility (which we may change in the next release\). Rychlik recommends a maximum value of 9 (Wojciech Rychlik, 'Selection of Primers for Polymerase Chain Reaction' in BA White, Ed., 'Methods in Molecular Biology, Vol. 15: PCR Protocols: Current Methods and Applications', 1993, pp 31-40, Humana Press, Totowa NJ\).)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# KM 8.11. 2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")


options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
primer3.path <- file.path(chipster.tools.path, "primer3" ,"src")

#check sequece file type
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, "sequence" )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "eprimer3")
emboss.parameters <- paste('sequence -auto -outfile primers.txt')
emboss.parameters <- paste(emboss.parameters, "-primer", primer)
emboss.parameters <- paste(emboss.parameters, "-task", task)
emboss.parameters <- paste(emboss.parameters, "-hybridprobe", hybridprobe)
emboss.parameters <- paste(emboss.parameters, "-numreturn", numreturn)
if (nchar(includedregion) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-includedregion", includedregion)
}
if (nchar(targetregion) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-targetregion", targetregion)
}
if (nchar(excludedregion) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-excludedregion", excludedregion)
}
if (nchar(forwardinput) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-forwardinput", forwardinput)
}
if (nchar(reverseinput) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-reverseinput", reverseinput)
}
emboss.parameters <- paste(emboss.parameters, "-gcclamp", gcclamp)
emboss.parameters <- paste(emboss.parameters, "-osize", osize)
emboss.parameters <- paste(emboss.parameters, "-minsize", minsize)
emboss.parameters <- paste(emboss.parameters, "-maxsize", maxsize)
emboss.parameters <- paste(emboss.parameters, "-opttm", opttm)
emboss.parameters <- paste(emboss.parameters, "-mintm", mintm)
emboss.parameters <- paste(emboss.parameters, "-maxtm", maxtm)
emboss.parameters <- paste(emboss.parameters, "-maxdifftm", maxdifftm)
emboss.parameters <- paste(emboss.parameters, "-ogcpercent", ogcpercent)
emboss.parameters <- paste(emboss.parameters, "-mingc", mingc)
emboss.parameters <- paste(emboss.parameters, "-maxgc", maxgc)
emboss.parameters <- paste(emboss.parameters, "-saltconc", saltconc)
emboss.parameters <- paste(emboss.parameters, "-dnaconc", dnaconc)
emboss.parameters <- paste(emboss.parameters, "-maxpolyx", maxpolyx)
emboss.parameters <- paste(emboss.parameters, "-psizeopt", psizeopt)
emboss.parameters <- paste(emboss.parameters, "-prange", prange)
emboss.parameters <- paste(emboss.parameters, "-ptmopt", ptmopt)
emboss.parameters <- paste(emboss.parameters, "-ptmmin", ptmmin)
emboss.parameters <- paste(emboss.parameters, "-ptmmax", ptmmax)
if (nchar(oexcludedregion) > 0 ) {
  emboss.parameters <- paste(emboss.parameters, "-oexcludedregion", oexcludedregion)
}
if (nchar(oligoinput) > 0 ) {
  emboss.parameters <- paste(emboss.parameters, "-oligoinput", oligoinput)
}
emboss.parameters <- paste(emboss.parameters, "-osizeopt", osizeopt)
emboss.parameters <- paste(emboss.parameters, "-ominsize", ominsize)
emboss.parameters <- paste(emboss.parameters, "-omaxsize", omaxsize)
emboss.parameters <- paste(emboss.parameters, "-otmopt", otmopt)
emboss.parameters <- paste(emboss.parameters, "-otmmin", otmmin)
emboss.parameters <- paste(emboss.parameters, "-otmmax", otmmax)
emboss.parameters <- paste(emboss.parameters, "-ogcopt", ogcopt)
emboss.parameters <- paste(emboss.parameters, "-ogcmin", ogcmin)
emboss.parameters <- paste(emboss.parameters, "-ogcmax", ogcmax)
emboss.parameters <- paste(emboss.parameters, "-osaltconc", osaltconc)
emboss.parameters <- paste(emboss.parameters, "-odnaconc", odnaconc)
emboss.parameters <- paste(emboss.parameters, "-oanyself", oanyself)
emboss.parameters <- paste(emboss.parameters, "-oendself", oendself)
emboss.parameters <- paste(emboss.parameters, "-opolyxmax", opolyxmax)
emboss.parameters <- paste(emboss.parameters, "-omishybmax", omishybmax)
emboss.parameters <- paste(emboss.parameters, "-explainflag", explainflag)
emboss.parameters <- paste(emboss.parameters, "-fileflag", fileflag)
emboss.parameters <- paste(emboss.parameters, "-firstbaseindex", firstbaseindex)
emboss.parameters <- paste(emboss.parameters, "-pickanyway", pickanyway)
emboss.parameters <- paste(emboss.parameters, "-maxmispriming", maxmispriming)
emboss.parameters <- paste(emboss.parameters, "-pairmaxmispriming", pairmaxmispriming)
emboss.parameters <- paste(emboss.parameters, "-numnsaccepted", numnsaccepted)
emboss.parameters <- paste(emboss.parameters, "-selfany", selfany)
emboss.parameters <- paste(emboss.parameters, "-selfend", selfend)
emboss.parameters <- paste(emboss.parameters, "-maxendstability", maxendstability)

eprimer.setup <- paste("export EMBOSS_PRIMER3_CORE=", primer3.path, "/primer3_core ;" ,sep="")
command.full <- paste(eprimer.setup, emboss.binary, emboss.parameters, ' >> eprimer3.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> eprimer3.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f eprimer3.log")
}
