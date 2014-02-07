#!/bin/bash -l
# K.M. 20.1. 2014
#
#
pblogfile=("./pb_blast_log")
pbtmproot=("/tmp")
start_time=$(date)
user=$(whoami)


if [[ $# -le 2 ]]
then
    echo "--------------------------------------------------------------------"
    echo "pb program is designed to run effectively BLAST searches"
    echo "for large amounts of query sequences."
    echo "pb accepts blast search commands and all their arguments."
    echo " "
    echo "Usage: pb blast_command e.g.:"
    echo "pb blastn -db embl -query myseqs.fasta -out results"
    echo " "
    echo "pb specific arguments:"
    echo "  -dbprot [seqfile]  Use the given protein sequce file as search database "
    echo "  -ensembl_dna [species_name]  Retrieve the genomic sequence of the given species from ensembl database and use it as the database"
    echo "  -ensembl_cdna [species_name]  Retrieve the coding sequences of the given species from ensembl database and use it as the database"
    echo "  -ensembl_prot [species_name]  Retrieve the peptide sequences of the given species from ensembl database and use it as the database"   
    echo "  -dbnuc [seqfile]  Use the given nucleotide sequce file as search database "
    echo "  -batch submit the job and write out instuctions how to follow the job and collect the results"
    echo "  -no_slurm Run the job in local machine and not through the batch job system"
    echo "--------------------------------------------------------------------"
    exit
fi

# Tarkista onko ohjelma sallittu Tarkista onko ohjelma sallittu
case "$1" in
       'blastp')
            blast_type=$1
        ;;
       'blastn')
            blast_type=$1
        ;;
       'blastx')
            blast_type=$1
        ;;
        'rpsblast')
            blast_type="blastp"
        ;;
       'psiblast')
            blast_type="blastp"
        ;;
        'rpsblastn')
            blast_type="blastx"
        ;;
       'tblastn')
            blast_type=$1
        ;;
       'tblastx')
            blast_type="blastn"
        ;;
       'deltablast')
            blast_type="blastp"
        ;;
        *)
            echo "   "
            echo "$argv[1] can not be used with pb"
            echo "pb can only be used to run commands:"
	        echo "  blastn"
		    echo "  blastp"
		        echo "  blastx"
            echo "  deltablast"
            echo "  psiblast"
            echo "  rpsblast"
            echo "  rpsblastn"
            echo "  tblastn"
            echo "  tblastx"
	        echo "  "
            echo "Usage: pb blast_command, e.g.:"
            echo "pb blastn -db embl -query myseqs.fasta -out results"
            exit 1
        ;;
esac



#module load emboss
command=(" ")
n_split=(500)


#Tarkista input, output ja tietokanta
inputflag=(0)
outputflag=(0)
databaseflag=(0)
ensemblflag=(0)
dbbuildflag=(0)
db_typy=(none)
pb_output_type=(0)
mode=("interactive")
debug=(0)
parse_seqids=("-parse_seqids")
blastdb_orig=($BLASTDB)
remote=(0)
num_threads=(1)

while [[ $# -ge 1 ]]
do
  case "$1" in
             '-query')
             # query file
                  inputseq=($2)
                  if [[ ! -e $inputseq ]]
                  then
                     echo ""
                     echo " Query sequence file: $inputseq does not exist"
                     echo ""
                     echo "-----------------------------------------------------------"
                     exit 1
                  fi
                  inputflag=(1)
			  inp_type=$(head -40 $inputseq | infoseq -nohead -only -type -filter | tail -1 | tr -d " " )
                  shift
                  shift
                ;;
               '-out')
               # result file
                  blast_output=($2)
                  outputflag=(1)
                  if [[ -e $blast_output ]]
                  then
                           echo " "
		           echo "Output file: $blast_output already exists."
                           echo ""
                           exit 1
                  fi
                  shift
                  shift
                ;;
                '-db')
                 # database
                  blast_database=($2)
                  databaseflag=(1)
                  if [[ -e $BLASTDB/$blast_database.nal ]]
                  then
		           db_type=(nuc)
				  fi
					  if [[ -e $BLASTDB/$blast_database.nsq ]]
                  then
		           db_type=(nuc)
				  fi
					  if [[ -e $BLASTDB/$blast_database.pal ]]
                  then
		           db_type=(prot)
				  fi
					  if [[ -e $BLASTDB/$blast_database.psq ]]
                  then
		           db_type=(prot)
				  fi
                  shift
                  shift
                ;;
                '-dbprot')
                 # database
                  blast_database=($2)
                  dbbuildflag=(1)
                  databaseflag=(1)
                  db_type=(prot)
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter -sformat asis | tail -1 | tr -d " \t")
                  if [[ $db_letter != "P" ]]
                  then
			      echo "ERROR: Protein database expected. "
                      echo "Your database sequence does not look like a protein sequence." ; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:" ; echo ""
                      head $blast_database
                      exit -1
                  fi
                  shift
                  shift
                ;;
			'-dbnuc')
                # database
                  blast_database=($2)
                  dbbuildflag=(2)
                  databaseflag=(1)
			  db_type=(nuc)
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter -sformat asis | tail -1 | tr -d " ")
                  if [[ $db_letter != "N" ]]
                  then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database
                      exit 1
                  fi
                  shift
                  shift
                ;;
                '-ensembl_dna')
                # database
                  blast_database=$(ensemblfetch.sh -type dna $2 | tee | tail -1 )
                  if [[ -e $blast_database ]]
                  then
                     echo "The genome data file $blast_database succesfully retrieved"
                  else
                     echo "The Ensembl data was not foud"
		     	       echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                  fi
                  dbbuildflag=(2)
                  databaseflag=(1)
                  ensemblflag=(1)
			  db_type=(nuc)
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 | tr -d " ")
                  if [[ $db_letter != "N" ]]
                  then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database
                      exit 1
                  fi
                  shift
                  shift
                ;;
                '-ensembl_cdna')
                # database
                   blast_database=$(ensemblfetch.sh -type cdna $2 | tail -1 )
                   if [[ -e $blast_database ]]
                   then
                     echo "The genome data file $blast_database succesfully retrieved"
                   else
                     echo "The Ensembl data was not foud"
		     	       echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                   fi
                   dbbuildflag=(2)
                   databaseflag=(1)
                   ensemblflag=(1)
			   db_type=(nuc)
                   db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 | tr -d " ")
                   if [[ $db_letter != "N" ]]
                   then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database
                      exit 1
                   fi
                   shift
                   shift
                 ;;

                '-ensembl_prot')
                 # database
                   blast_database=$(ensemblfetch.sh -type pep $2 | tail -1 )
                   if [[ -e $blast_database ]]
                   then
                    echo "The protein data file $blast_database succesfully retrieved"
                   else
                     echo "The Ensembl data was not foud"
		     	       echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                   fi
                   dbbuildflag=(1)
                   databaseflag=(1)
                   ensemblflag=(1)
                   db_type=(prot)
                   db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 |  tr -d " " )
                   if [[ $db_letter != "P" ]]
                   then
		      echo "ERROR: Protein database expected. "
                      echo "Your database sequence does not look like a protein sequence." ; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:" ; echo ""
                      head $blast_database
                      exit
                   fi
                   shift
                   shift
                ;;
                 '-entrez_query')
                 command="$command $1 \"$2\""
                 shift
                 shift
                 ;;
		 '-outfmt')
                 f_len=$(echo $2 | wc -c)
                 #If the argument is longer than 4 then we assume that this is fromat 6,7 or 10 with definitions
                 if [[ $f_len -le 4 ]]
                 then
                   if [[ $2 -eq 12 || $2 -eq 13 ]]
                   then
                     fromatstring=('"6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore"')
                     command="$command $1 $fromatstring"
                     pb_output_type=($2)
                   fi
                   if [[ $2 -eq 14 ]]
                   then
                     fromatstring=('"6 sacc sstart send sseq"')
                     command="$command $1 $fromatstring "
                     pb_output_type=($2)
                   fi
                   if [[ $2 -eq 15 ]]
                   then
                      fromatstring=("10")
                      command="$command $1 $fromatstring "
                      pb_output_type=($2)
                    fi
                    if [[ $2 -eq 16 ]]
                    then
                      fromatstring=("6")
                      command="$command $1 $fromatstring  "
                      pb_output_type=($2)
                    fi
                    if [[ $2 -lt 12 ]]
                    then
                       command="$command $1 \"$2\""
                       pb_output_type=($2) 
                       echo $command
                    fi 
                 else
                     #parsing format types 6,7,10                    
                     formatarray=($2)
                     pb_output_type=${formatarray[0]}
                     format_n_items=${#formatarray[*]}
                     command="$command $1 \"$pb_output_type"
                     for ((item=1; item<$format_n_items; item++))
                     do
                       command="$command ${formatarray[$item]}"
                     done
                     command="$command\" "
                 fi
                 shift
                 shift
                ;;
                '-num_threads')
                    num_threads=($2)
                    shift
                    shift
                ;;

                '-chipster_path')
	           chipster_path="$2"
                   tools_path="$chipster_path""/tools"
                   comp_path="$chipster_path""/comp"
                   export PATH=${PATH}:$comp_path/modules/admin/shell/:$tools_path/emboss/bin/:$tools_path/blast/bin/
                   shift
                   shift
                ;;


                '-help')
			  $@
                  echo " *** pb spesific options:"
                  echo "  -dbprot  Use the given protein sequce file as search database "
                  echo "  -dbnuc  Use the given nucleotide sequce file as search database "
                  echo "  -ensembl_dna [species_name]  Retrieve the genomic sequence of the given species from ensembl database and use it as the database"
                  echo "  -ensembl_cdna [species_name]  Retrieve the coding sequences of the given species from ensembl database and use it as the database"
                  echo "  -ensembl_prot [species_name]  Retrieve the peptide sequences of the given species from ensembl database and use it as the database"
                  echo "  -batch submit the job and write out instuctions how to follow the job and collect the results"
                  echo "  -no_seqid_parsing  With this option the database indexing for a sequence set defined with -dbnuc or -dbprot is done without -parse_seqids option"                            echo " "
                  echo "Extra modes for -outfmt option"
                  echo "  -outfmt 12  Print out a list of uniq hit sequece identifiers" 
                  echo "  -outfmt 13 print the hit sequences in fasta format"
                  echo "  -outfmt 14  Print out the matching regions of hit sequences in fasta format" 
                  echo " "
                  exit
			;;
                '-nsplit')
                  n_split=($2)
                  shift
                  shift
                ;;
                '-batch')
                  mode=("batch")
                  shift
                ;;
                '-no_slurm')
                  mode=("no_slurm")
                  shift
                ;;
                '-debug')
                  debug=(1)
                  shift
                ;;
                '-no_seqid_parsing')
                  parse_seqids=(" ")
                shift
                ;;
                '-remote')
                  remote=(1)
                  command="$command $1"
                shift
                ;;
                '-archive')
                 inputseq=($2)
                 blast_database=($2)
                  inputflag=1
                  inp_type=$(file $inputseq | awk -F : '{print $2}')                              db_type=(archive)
                shift
                shift
                ;;
                *)
            	command="$command $1"
                shift                       # No more switches
                ;;
    esac
done



if [[ $inputflag == 0 ]]
then
   echo "   "
   echo "Inputfile not defined"
   echo "You must define the file containing input sequences with option: -query file_name"
   echo "  "
   exit 1
fi

if [[ $databaseflag == 0 || $db_type == "none" ]]
then
   echo "   "
   echo "Database not defined"
   echo "You must define the query database with option:"
   echo "          -db database_name"
   echo "or feed in fasta formatted sequence file to be used as database with options:" 
   echo "          -dbprot protein_sequence_file"
   echo "or "
   echo "          -dbnuc nucleotide_sequece_file"
   echo "or you can use options -ensembl_dna, -ensembl_cdna, -ensembl_prot to retrieve "
   echo "a target genome from the Ensembl databases"
   exit 1
fi

if [[ $blast_type == "blastn" ]] 
then
     if [[ $db_type == "prot" ]]
     then
     echo "    "
     echo "Wrong database type."
     echo "blastn can only be used for searches against nucleotide databases."
     echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]]
     then
     echo "    "
     echo "Wrong query sequence type."
     echo "blastn can only be used for searches with nucleotide sequences."
     echo "    "
        exit 1
     fi
fi

if [[ $blast_type == "blastp"  ]] 
then
     if [[ $db_type == "nuc" ]]
     then
     echo "    "
     echo "Wrong database type."
     echo "blastp can only be used for searches against protein databases."
     echo "    "
        exit 1
     fi
     if [[ $inp_type == "N" ]] 
     then
     echo "    "
     echo "Wrong query sequence type."
     echo "blastp can only be used for searches with protein sequences."
     echo "    "
        exit 1
     fi
fi

if [[ $blast_type == "blastx" ]] 
then
     if [[ $db_type == "nuc" ]]
     then
     echo "    "
     echo "Wrong database type."
     echo "blastx can only be used for searches against protein databases."
     echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]] 
     then
     echo "    "
     echo "Wrong query sequence type."
     echo "blastx can only be used for searches with nucleotide sequences."
     echo "    "
        exit 1
     fi
fi
if [[ $blast_type == "tblastn" ]] 
then
     if [[ $db_type == "prot" ]] 
     then
     echo "    "
     echo "Wrong database type."
     echo "tblastn can only be used for searches against nucleotide databases."
     echo "    "
        exit 1
    fi
     if [[ $inp_type == "N" ]] 
     then
     echo "    "
     echo "Wrong query sequence type."
     echo "tblastn can only be used for searches with protein sequences."
     echo "    "
        exit
     fi
fi

if [[ $blast_type == "tblastx" ]] 
then
     if [[ $db_type == "prot" ]]
     then
     echo "    "
     echo "Wrong database type."
     echo "tblastx can only be used for searches against nucleotide databases."
     echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]] 
     then
     echo "    "
     echo "Wrong query sequence type."
     echo "tblastx can only be used for searches with nucleotide sequences."
     echo "    "
        exit 1
     fi
fi

if [[ $outputflag == 0 ]]
then
     blast_output=(blast_results)
fi



echo "Create a temporary directory"
echo "mkdir $pbtmproot/pb_$$_tmpdir"

mkdir $pbtmproot/pb_$$_tmpdir


#Move and check the inputseq
cp $inputseq $pbtmproot/pb_$$_tmpdir/
inputname=$(ls  $pbtmproot/pb_$$_tmpdir/)
inputseq=($pbtmproot/pb_$$_tmpdir/$inputname)


#In Batch mode, create a temporaray version of database to avoid
#problems caused by the database update
if [[ $mode == "batch" && $remote == 0 ]]
then
  if [[ $dbbuildflag -lt 1 ]]
  then
      if [[ $outputflag -ne 0 ]]
      then
        echo "  "
        echo "Making temporary copy of $blast_database."
        echo " "
      fi
      blastdb_orig=($BLASTDB)

      # check if the database refres to a .pal or .nal file
      if [[ -e  "$BLASTDB"/"$blast_database".pal ]]
      then
       cp  "$BLASTDB"/"$blast_database".pal $pbtmproot/pb_$$_tmpdir/
          for subdb in $(grep DBLIST "$BLASTDB"/"$blast_database".pal | sed s/"DBLIST "/""/ )
          do
            location=$(pwd)
            cd $BLASTDB
            cp "$subdb".*  $pbtmproot/pb_$$_tmpdir/
            cd $location
          done
      else
         if [[ -e  "$BLASTDB"/"$blast_database".nal ]]
         then
           cp "$BLASTDB"/"$blast_database".nal $pbtmproot/pb_$$_tmpdir/
           for subdb in $(grep DBLIST "$BLASTDB"/"$blast_database".nal | sed s/"DBLIST "/""/)
           do
             location=$(pwd)
             cd $BLASTDB
             cp "$subdb".*  $pbtmproot/pb_$$_tmpdir/
             cd $location
           done
           else
           cp "$BLASTDB"/"$blast_database".*  $pbtmproot/pb_$$_tmpdir/
        fi
      fi
      export BLASTDB=$pbtmproot/pb_$$_tmpdir 
      #ls $BLASTDB
  fi
fi

#Build BLAST-index if needed

#First check and fix possible duplicate database entries
if [[ $dbbuildflag -gt 0 ]]
then
  grep ">" $blast_database | cut -f1 -d " " | sort | uniq -c | awk '{if ($1 > 1) print $0}'  >  $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names
  nonuniq=$(cat  $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names | wc -l )
  if [[ $nonuniq -gt 0 ]]
  then
    echo "Following sequence names existed in the database sequence file more than once"
    cat $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names | awk '{print $2 " ocurrences: " $1}' 
    echo "Fixing names"
    cp $blast_database $blast_database"_tmp_"$$
    for non_uniq_name in $(awk '{print $2}' $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names)
    do
        echo "fixing $non_uniq_name"
        awk 'BEGIN{n=1}{ if ( $1 == "'"$non_uniq_name"'" ) print $1"_"n" "$2" "$3" "$4" "$5" "$6} {if ( $1 == "'"$non_uniq_name"'" ) n = n +1}{ if ( $1 != "'"$non_uniq_name"'" ) print $0} ' $blast_database"_tmp_"$$ >>!  "$blast_database"_"$$"_fixing
        rm -f $blast_database"_tmp_"$$
        mv "$blast_database"_"$$"_fixing $blast_database"_tmp_"$$
     done
     blast_database=("$blast_database"_tmp_"$$")
  fi
fi


if [[ $dbbuildflag -gt 0 ]]
then

  if [[ $dbbuildflag -eq 1 ]]
  then
    dbtype=(prot)
  elif [[ $dbbuildflag -eq 2 ]]
  then
    dbtype=(nucl)
  fi

  if [[ $outputflag -ne 0 ]]
  then
    echo "  "
    echo "Building BLAST indexes for sequences in $blast_database."
    echo " "
  fi

  num_db_seq=$(grep -c ">" $blast_database)


  #salloc -n 1 -t 8:00:00
  #srun seqret $blast_database $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta
  #seqret $blast_database $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta
  cp $blast_database $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta

  makeblastdb -dbtype $dbtype -in $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta $parse_seqids -out $pbtmproot/pb_$$_tmpdir/pb_$$_blast_database.tmp  > $pbtmproot/pb_$$_tmpdir/makeblastdb_log

  blastdb_num=$(grep "Adding sequences from" $pbtmproot/pb_$$_tmpdir/makeblastdb_log | awk '{print $6}')
  if [[ $num_db_seq -ne $blastdb_num ]]
  then
    cat $pbtmproot/pb_$$_tmpdir/makeblastdb_log
    exit 1
  elif [[ $outputflag -ne 0 ]]
  then
     echo "$num_db_seq sequences indexed from $blast_database"
  fi

  blastdb_orig=($BLASTDB)
  export BLASTDB=$pbtmproot/pb_"$$"_tmpdir
  blast_database_name=($blast_database)
  blast_database=(pb_"$$"_blast_database.tmp)
fi

# Running the job
if [[ $outputflag -ne 0 ]]
then
   echo "  "
   echo "Running:"
   if [[ $blast_type == "blast_formatter" ]]
   then
      echo $command -archive $inputseq -out $blast_output
   else
      echo "$command -query $inputseq -db $blast_database -out $blast_output"
   fi
   echo " "
fi


#Syötedatan käsittely tehdään eri tavalla jos kyse on archive-tiedostosta
if [[ $blast_type != "blast_formatter" ]]
then
   if [[ $inp_format != "fasta" ]]
   then
      echo "Converting $inp_format formated query file to fasta format with seqret"
      seqret $inputseq "$inputseq".fasta >& /dev/null
      inputseq=(""$inputseq".fasta")
   fi
   infoseq $inputseq -only -usa -nohead -outfile ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list >& /dev/null
else
   # datan käsittely blast_formatterissa
   location=$(pwd)
   cd $pbtmproot/pb_"$$"_tmpdir/
   tar xf $inputseq
   num=$(ls *blast_database* | wc -l )
   if [[ $num > 1 ]]
   then
      export BLASTDB=${pbtmproot}/pb_"$$"_tmpdir/
   fi
   cd $location
   ls ${pbtmproot}/pb_"$$"_tmpdir/*result > ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list
   cp ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_chunk_list
fi


r=$(cat ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list | wc -l)

if [[ $r -eq 0 ]]
then
   echo "  "
   echo "Error: Input sequence $inputseq not found"
   exit 1
fi



koko=$(expr $r / $n_split)
if [[ $blast_type == "blast_formatter" ]]
then
     # Kaikki OK
   echo "blast_formatter"
else
   # No splitting in no_slurm mode
   if [[ $mode == "no_slurm" ]]
   then
     cp $inputseq  ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_00001.fasta
   else
     echo "split"
     if [[ $koko -gt 5 ]]
     then
         split_fasta.pl -p ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_ -s .fasta -o 1 $koko <  $inputseq  >& /dev/null
     else
         split_fasta.pl -p ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_ -s .fasta -o 1 5 <  $inputseq >& /dev/null
     fi
   fi
   ls ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_*.fasta > $pbtmproot/pb_$$_tmpdir/pb_$$_chunk_list
fi

n_split=$(cat ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_chunk_list | wc -l )

if [[ $outputflag -ne 0 ]]
then
    echo "Total amount of query sequences: $r"
    echo "The job is split into $n_split pieces"
fi

location=$(pwd)
cd ${pbtmproot}/pb_"$$"_tmpdir/


#######
#Kirjoitetaan era-ajotiedosto
######

#check the database size

if [[ -e $BLASTDB/$blast_database.pal ]]
then
    koko=(0)
    for osa in $(grep DBLIST $BLASTDB/$blast_database.pal | cut -c7-900 | sed s/"\/qfs2\/biodb\/blast\/"//g )
    do
      osa=$( echo $osa | awk -F "/" '{ print $NF}' )
      lisa=$( ls -l ${BLASTDB}/$osa* | awk '{ a = (a + $5)}{print a}' | tail -1 )
      (( koko = koko + lisa ))
    done
elif [[ -e $BLASTDB/$blast_database.nal ]]
then
    koko=(0)
    for osa in $(grep DBLIST $BLASTDB/$blast_database.nal | cut -c7-900 | sed s/"\/qfs2\/biodb\/blast\/"//g ) 
    do
      osa=$(echo $osa | awk -F "/" '{ print $NF}' )
      lisa=$(ls -l $BLASTDB/$osa* | awk '{ a = (a + $5)}{print a}' | tail -1 )
      koko=$(expr $koko + $lisa)
    done
elif [[ $remote == 0 ]]
then
    koko=$(ls -l $BLASTDB/$blast_database* | awk '{ a = (a + $5)}{print a}' | tail -1)
else
    koko=6990510
fi

# 2 * 1048576 / 3 = 699051
(( mem_request = koko / 699051 ))

if [[ $mem_request -lt 2000 ]]
then
  mem_request=2000
fi

if [[ $mem_request -gt 64000 ]]
then
  mem_request=64000
fi

if [[ $outputflag != 0 ]]
then
   echo "requesting memory $mem_request MB"
fi

#Set the number of cores:
if [[ $blast_type == "blast_formatter" ]]
then
  num_of_cores_to_use=(1)
else
  num_of_cores_to_use=($num_threads)
fi

(( mem_request = mem_request / num_of_cores_to_use ))

#
#Create the batch job script
#
echo '#!/bin/bash' > $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -J pb_blast' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -o pb_blast_out_%A_%a.txt' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -e pb_blast_err_%A_%a.txt' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -t 24:00:00' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --mem-per-cpu='"$mem_request" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --array=1-'"$n_split" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -n 1' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --cpus-per-task='"$num_of_cores_to_use" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo ""  >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
#echo "cd $pbtmproot/pb_"$$"_tmpdir/jobs/job_\$SLURM_ARRAY_TASK_ID" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo "export BLASTDB=$BLASTDB" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

if [[ $mode == "no_slurm" ]]
then
echo "SLURM_ARRAY_TASK_ID=1" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
fi


echo 'nimi=(`sed -n "$SLURM_ARRAY_TASK_ID"p $pbtmproot/pb_'$$'_tmpdir/pb_'$$'_chunk_list `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'start_time=(`date +%s`)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
if [[ $blast_type == "blast_formatter" ]]
then
    echo $command '-archive $nimi -out $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

else
    if [[ $mode == "no_slurm" ]]
    then
        echo $command '-query $nimi -db' $blast_database '-out $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    else
        echo $command '-num_threads $SLURM_CPUS_PER_TASK -query $nimi -db' $blast_database '-out $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    fi
fi

#  #add missing newline character to XML formatted outputfiles
#  if ( $pb_output_type == 5) then 
#      echo 'echo "" >> $nimi.result' >>  $pbtmproot/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#  endif

 # ID listan luonti tuloksista
if [[ $pb_output_type -eq 12 || $pb_output_type -eq 13 ]] 
then    
    echo 'numhits=$(cat $nimi.result | wc -l )' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'if [[ $numhits -gt 0 ]] ; then ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  awk '\''{print $2}'\'' $nimi.result | sort | uniq  > $nimi.result.list' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  rm -f $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  mv $nimi.result.list $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'fi' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

fi 


#if ( $pb_output_type == 12) then
#    echo 'rm -f $nimi.result' >> $pbtmproot/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#    echo 'mv $nimi.result.list $nimi.result' >> $pbtmproot/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#endif
  
echo 'end_time=(`date +%s`)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '((duration = end_time - start_time))' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'echo '$jobindex' Duration: $duration'  >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh


if [[ $pb_output_type -eq 14 ]] 
then 
  echo 'numhits=(`cat $nimi.result | wc -l `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo 'if [[ $numhits -gt 0 ]]; then  ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  sort $nimi.result | uniq > $nimi.result.uniq' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  rm -f $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  mv $nimi.result.uniq $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
  echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

if [[ $pb_output_type -eq 15 ]]
then
   echo 'numhits=(`cat $nimi.result | wc -l `)' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'if [[ $numhits -gt 0 ]] ; then  ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  awk -F "," '\''{print $2}'\'' $nimi.result > $nimi.list' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  blastdbcmd -db '$blast_database' -entry_batch $nimi.list -outfmt %t -out $nimi.desc ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  perl /appl/bio/blast/ncbi-blast-2.2.28+/bin/yhd.pl $nimi.result $nimi.desc COMMA > $nimi.ext ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  rm -f $nimi.result ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  mv $nimi.ext $nimi.result' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

if [[ $pb_output_type -eq 16 ]]
then
   echo 'numhits=(`cat $nimi.result | wc -l `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'fi [[ $numhits -gt 0 ]]; then  ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  awk '\''{print $2}'\'' $nimi.result > $nimi.list' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  blastdbcmd -db '$blast_database' -entry_batch $nimi.list -outfmt %t -out $nimi.desc ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  perl  /appl/bio/blast/ncbi-blast-2.2.28+/bin/yhd.pl $nimi.result $nimi.desc TAB > $nimi.ext ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  rm -f $nimi.result ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  mv $nimi.ext $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

#############
#Tyon lahetys
##############
cd  $pbtmproot/pb_"$$"_tmpdir/

if [[ $mode == "no_slurm" ]] 
then
  source pb_batch.sh > pb_blast_out 2> pb_blast_err
else 
 jobid=$(sbatch pb_batch.sh | grep "Submitted" | awk '{print $NF}')
 echo $jobid

 n_done=(0)
 n_failed=(0)
 n_finished=(0)
 n_running=(0)
 n_waiting=(0)
 i=(0)

 while [[ $n_finished -lt $n_split ]]
 do
    prew_running=($n_running)
    prew_waiting=($n_waiting)
    prew_done=($n_done)
    n_running=$(scontrol show jobs $jobid | grep JobState | grep RUNNING | wc -l )
    n_done=$(scontrol show jobs $jobid | grep JobState | grep COMPLETED | wc -l )
    n_waiting=$(scontrol show jobs $jobid | grep JobState | grep PENDING | wc -l )
    n_failed=$(scontrol show jobs $jobid | grep JobState | grep FAILED | wc -l )
    (( n_finished =  n_split - n_waiting - n_running ))
    change=(0)
      if [[ $prew_running != $n_running ]]
      then 
        change=(1)
      fi
      if [[ $prew_waiting != $n_waiting ]] 
      then 
        change=(1)
      fi
      if [[ $change == 1 ]] 
      then  
          if [[ $i -gt 0 ]] 
          then
            echo ""
          fi  
          aika=$(date +%H:%M:%S)
          printf "$aika Jobs: done $n_finished, running $n_running, waiting $n_waiting, failed $n_failed\n"
      else 
        if [[ $i -lt 80 ]] 
        then
         printf "%s" .
         (( i = i + 1 ))
       else 
         printf "%s\n" .
         i=(0)
        fi 
     fi 
   sleep 15
 done
fi

cd $location

#Kootaan lopputulokset

#Jos käytetään arkistoformaattia, tulostiedostoja ei suoraan voida yhdistää.
if [[ $pb_output_type == 11 ]] 
then
    cd $pbtmproot/pb_$$_tmpdir/
    tar cf $blast_output.tar pb_chunk_*.result
    if [[ $dbbuildflag -eq 2 || $dbbuildflag -eq 1 ]] 
    then
       tar -rf $blast_output.tar  pb*blast_database.tmp*
    fi
    mv $blast_output.tar $location/
    cd $location
    blast_output=("$blast_output.tar")
else 
    if [[ $blast_type == "blast_formatter" ]]
    then 
       cat $pbtmproot/pb_$$_tmpdir/pb_chunk_*.result.result > $blast_output
    else
       cat $pbtmproot/pb_$$_tmpdir/pb_chunk_*.result > $blast_output
    fi
fi


# Nimilista: poistetaan päällekäisyydet
if [[ $pb_output_type -eq 12 || $pb_output_type -eq 13 ]] 
then
  sort $blast_output | uniq > "$blast_output"_$$.tmp
  rm -f $blast_output
  mv "$blast_output"_$$.tmp $blast_output
    numhits=$(cat $blast_output | wc -l )
fi

# Haetaan sekvenssit
if [[ $pb_output_type -eq 13 ]] 
then
  if [[ $numhits -gt 0 ]] 
  then
    blastdbcmd -db $blast_database -entry_batch $blast_output -out $blast_output.seq
    rm -f $blast_output
    seqret $blast_output.seq $blast_output
    rm -f $blast_output.seq
  fi
fi

# Muokataan taulukko fasta-tiedostoksi
if [[ $pb_output_type -eq 14 ]] 
then 
  numhits=$(cat $blast_output | wc -l )
  if [[ $numhits -gt 0 ]] 
  then
     # aukkojen erilaisuuden takia listassa saatta olla vielä päällekäisyyttä
          sort $blast_output | awk '{print $4"\t"$1"\t"$2"\t"$3}' | uniq -f 1 | awk -F "\t" '{print ">"$2"_"$3"-"$4 }{print $1}' | degapseq -filter > "$blast_output"_"$$".seq
     rm -f $blast_output
     mv "$blast_output"_"$$".seq $blast_output
  fi
fi

# Muokataan taulukkoa ja lisätään otsikkorivi
if [[ $pb_output_type -eq 16 ]] 
then 
  numhits=$(cat $blast_output | wc -l )
  if [[ $numhits -gt 0 ]] 
  then 
     echo "query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tq start\tq. end\ts. start\ts. end\te-value\tbit score\ts.description" >"$blast_output"_$$.table
     cat $blast_output >> "$blast_output"_$$.table
     rm -f $blast_output
     mv "$blast_output"_$$.table $blast_output
  fi 
fi

if [[ $outputflag -eq 0 ]] 
then  
  cat $blast_output
  rm -f $blast_output
fi

if [[ $outputflag -eq 1 ]] 
then  
   aika=$(date +%H:%M:%S)
   printf "$aika Job completed. Results are written to file: $blast_output\n"
fi


if [[ $dbbuildflag -ne 0 ]] 
then
   export BLASTDB=$blastdb_orig
fi

cat $pbtmproot/pb_$$_tmpdir/pb_blast_err* > $pbtmproot/pb_$$_tmpdir/pb_errors
grep "Duration:" $pbtmproot/pb_$$_tmpdir/pb_blast_out* > subjobtimes_$$

if [[ -e  $pbtmproot/pb_$$_tmpdir/pb_errors ]] 
then 
  cat $pbtmproot/pb_$$_tmpdir/pb_errors
fi

 
if [[ $ensemblflag == 1 ]] 
then 
   rm -rf $blast_database_name
fi

rm -rf $pbtmproot/pb_$$_tmpdir 
#rm -f pb_blast_job_tmp_$$* 
end_time=$(date)

if [[ $ensemblflag == 1 ]] 
then
    blast_database=($blast_database_name)
    rm -rf $blast_database_name
fi

#echo $user $start_time $end_time $blast_database $r >> $pblogfile

exit 0
  