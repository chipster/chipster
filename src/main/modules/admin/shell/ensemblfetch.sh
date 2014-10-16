#!/bin/bash 
# A script to fetch datasets from the ensembl website
# 23.7. 2010 KM
# Updated 4.11. 2011 KM

TMPDIR="/tmp"

if [[ "$1" == "" ]]
 then
   echo give a species name
   exit
fi

seqtype=("dna")
mode=("single")
outputmode=("single")
ensembl_version=("current_")
ensembl_genomes_version=("current")
ensembl_urls_file="$TMPDIR/ensembl_urls_${ensembl_genomes_version}"
ensembl_urls_mysql_file="$TMPDIR/ensembl_urls_mysql_${ensembl_genomes_version}"

update_ensembl_url_file ()
{
  url_file=$1
  type=$2 # fasta or mysql
  
  if [[ ! -e $url_file ]]
  then
    # create a list of ftp directories
    echo ftp://ftp.ensembl.org/pub/$ensembl_version$type/ > ensembl_list

    # mysql data for bacterisa is in different format
    if [[ "$type" != "mysql" ]]
    then
      for ((number=1; number<=27; number++))
      do
        echo ftp://ftp.ensemblgenomes.org/pub/bacteria/${ensembl_genomes_version}/$type/bacteria_"$number"_collection/ >> ensembl_list
      done
    fi

    echo ftp://ftp.ensemblgenomes.org/pub/fungi/${ensembl_genomes_version}/$type/ >> ensembl_list
    echo ftp://ftp.ensemblgenomes.org/pub/metazoa/${ensembl_genomes_version}/$type/  >> ensembl_list
    echo ftp://ftp.ensemblgenomes.org/pub/plants/${ensembl_genomes_version}/$type/ >> ensembl_list
    echo ftp://ftp.ensemblgenomes.org/pub/protists/${ensembl_genomes_version}/$type/  >> ensembl_list

    # create a list of species in all directories
    rm -f ensembl_urls
    while read url
    do
      # there may be a proxy (squid) set in environment variable to make file downloads faster,
      # but it must be disabled for directory listing, because squid converts directory listing to html
      curl --proxy "" --silent --list-only $url > directories
      while read directory
      do
        echo "$url$directory/" >> ensembl_urls
      done < directories
    done < ensembl_list

    # remove duplicates, because some species are on multiple sites
    cat ensembl_urls \
      | grep -v "ensembl.org/pub/${ensembl_version}fasta/caenorhabditis_elegans" \
      | grep -v "ensembl.org/pub/${ensembl_version}fasta/saccharomyces_cerevisiae" \
      | grep -v "ensemblgenomes.org/pub/metazoa/${ensembl_genomes_version}/fasta/drosophila_melanogaster" \
      > filtered_urls

    cat filtered_urls > $url_file
    if [[ "$type" == "mysql" ]]
    then
      cat filtered_urls | grep "_core_" > $url_file
    else
      cat filtered_urls > $url_file
    fi
  fi
}

update_ensembl_urls () 
{
  echo "Getting the name list"
  update_ensembl_url_file $ensembl_urls_file "fasta"
  update_ensembl_url_file $ensembl_urls_mysql_file "mysql"
}

clean_up ()
{
  rm -rf  tmp_$$/
}

while [[ $# -ge 1 ]]
do
  case "$1" in
             #data type
             '-type')
                seqtype="$2"
                shift
                shift
              ;;
              #
              '-names')
                  echo "Retrieving the list of available species names:"

                   rm -f $ensembl_urls_file
                   rm -f $ensembl_urls_mysql_file

                   mkdir tmp_$$
                   cd tmp_$$
                   update_ensembl_urls
                   cd ..

                   cat $ensembl_urls_file |\
                     sed s/"bacillus_collection\/b_"/"bacillus_"/g | \
                     sed s/"borrelia_collection\/b_"/"borrelia_"/g | \
                     sed s/"buchnera_collection\/b_"/"buchnera_"/g | \
                     sed s/"escherichia_shigella_collection\/e_"/"escherichia_"/g | \
                     sed s/"escherichia_shigella_collection\/s_"/"shigella_"/g | \
                     sed s/"mycobacterium_collection\/m_"/"mycobacterium_"/g | \
                     sed s/"neisseria_collection\/n_"/"neisseria_"/g | \
                     sed s/"pyrococcus_collection\/p_"/"pyrococcus_"/g |\
                     sed s/"streptococcus_collection\/s_"/"streptococcus_"/g |\
                     sed s/"staphylococcus_collection\/s_"/"staphylococcus_"/g |\
                   awk -F "/" '{print $(NF-1)}' | sort
                   clean_up
                   exit
               ;;
              '-list')
                  namelist=($2)
                  spec="species in file $namelist"
                  if ![[ -e "$namelist" ]]
                  then
                     echo "Species name list file: $namelist not found!"
                     exit
                  fi
                  shift
                  shift
                  mode=("list")
              ;;
              '-version')
                  ensembl_version=("release-$2/")
                  ensembl_genomes_version=("release-$2")  
                  ensembl_urls_file="$TMPDIR/ensembl_urls_${ensembl_genomes_version}"
                  ensembl_urls_mysql_file="$TMPDIR/ensembl_urls_mysql_${ensembl_genomes_version}"
		  echo "Using Ensembl version $ensembl_genomes_version"
                  shift
                  shift
              ;;
              '-out')
               outfile=($2)
                  if [[ -e "$outfile" ]]
                  then
                     echo "Result file: $outfile already exists!"
                     exit
                  fi
                  shift
                  shift
                  outputmode=("file")
              ;;
              '-tmp')
                 TMPDIR=($2)
                 ensembl_urls_file="$TMPDIR/ensembl_urls_${ensembl_genomes_version}"
                 ensembl_urls_mysql_file="$TMPDIR/ensembl_urls_mysql_${ensembl_genomes_version}"
                 shift
                 shift
              ;;
              '-help')
                 echo " ----------------------------------------------------------------------------------"
                 echo " ensemblfetch retrieves genomic, cDNA or peptide sequences of a given species from "
                 echo " the ensembl or ensembl genomes ftp site."
                 echo " "
                 echo " Usage:"
                 echo "   ensemblfetch [-options] species_name"
                 echo ""
                 echo " ensemblfetch options:"
                 echo "  -names   List the available species names "
                 echo "  -type    Select the data type to retrieve (Default: dna) "
                 echo "           Available data types are: dna, dna_rm, cdna, cdna_abinitio, pep, pep_abinitio"
                 echo "  -list    List of species to retrieve"
                 echo "  -version Download from specified ensembl release instead of latest"
                 echo "  -out     Output file name" 
                 echo "  -help    print this help"
                 echo " "
                 exit
              ;;
              *)
                spec=($1)
                shift                       # No more switches
              ;;
    esac
done
case "$seqtype" in
    "dna")
       echo "Retrieving genomic DNA for $spec"
     ;;
     "dna_rm")
       echo "Retrieving repeat masked genomic DNA for $spec"
     ;;
     "cdna")
       echo "Retrieving sequences for all transcript (known, novel and pseudo) for $spec" 
     ;;
     "cdna_abinitio")
       echo "Retrieving sequences for ab-initio prediction transcripts for $spec"
     ;;
     "pep")
       echo "Retrieving sequences for all all known and novel peptides for $spec"
     ;;
     "pep_abinitio")
       echo "Retrieving sequences for all abinitio predicted peptides for $spec"
     ;;
    "gtf")
       echo "Retrieving gtf file for $spec"
     ;;
    "mysql")
       echo "Retrieving mysql files for $spec"
     ;;
     *)
         echo "Unknown data type"
         echo "Please use one of the following types:"
         echo "  dna"
         echo "  dna_rm"
         echo "  cdna"
         echo "  cdna_abinitio"
         echo "  pep"
         echo "  pep_abinitio"
         echo "  gtf"
         echo "  mysql"
         exit
     ;;
esac


mkdir tmp_$$
cd tmp_$$

update_ensembl_urls

if [[ $seqtype == "dna" ]]
then
   awk '{print $1"dna/*.dna.toplevel.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "dna_rm" ]]
then
   awk '{print $1"dna/*.dna_rm.toplevel.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "cdna" ]]
then
   awk '{print $1"cdna/*.cdna.all.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "cdna_abinitio" ]]
then
   awk '{print $1"cdna/*.cdna.abinitio.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "pep" ]] 
then
   awk '{print $1"pep/*.pep.all.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "pep_abinitio" ]]
then
   awk '{print $1"pep/*.pep.abinitio.fa.gz"}' $ensembl_urls_file > ensembl_species.txt
fi

if [[ $seqtype == "gtf" ]]
then
   awk '{print $1"*gtf.gz"}' $ensembl_urls_file | sed s/"fasta"/"gtf"/g > ensembl_species.txt
fi

if [[ $seqtype == "mysql" ]]
then
   cat $ensembl_urls_mysql_file > ensembl_species.txt
fi
 
cd ..

if [[ $mode == "single" ]]
then 
   echo "$spec" | tr "[:upper:]" "[:lower:]" | sed s/" "/"_"/g  > tmp_$$/namelist
else
   cat $namelist | tr "[:upper:]" "[:lower:]" | sed s/" "/"_"/g > tmp_$$/namelist
fi

for species in $(cat tmp_$$/namelist)
do
  if [[ $seqtype == "mysql" ]]
  then
     numhits=$(grep -i "/"$species"_core" tmp_$$/ensembl_species.txt| wc -l)
  else
     numhits=$(grep -i "/$species/" tmp_$$/ensembl_species.txt| wc -l)
  fi
  if [[ $numhits -eq 0 ]]
  then
    echo "--------------------------------------------------------------------------------"
    echo "The species name: $species was not found from "
    echo "the ensembl and ensembl genomes databases"
    echo "--------------------------------------------------------------------------------"
  fi

  if [[ $numhits -eq 1 ]]
  then
    if [[ $seqtype == "mysql" ]]
    then
      url=$(grep -i "/$species" tmp_$$/ensembl_species.txt)
      filename=$(grep -i "/$species" tmp_$$/ensembl_species.txt | awk -F "/" '{print $(NF-1)}' )
      mkdir "$filename"_mysql
      cd "$filename"_mysql
      wget -o log "$url""coord_system.txt.gz" >> /dev/null
      wget -o log "$url""seq_region.txt.gz" >> /dev/null
      wget -o log "$url""karyotype.txt.gz" >> /dev/null
      wget -o log "$url""repeat_feature.txt.gz" >> /dev/null
      cd ..
      tar cvf "$filename"_mysql.tar "$filename"_mysql
      rm -rf "$filename"_mysql
      echo "The results have been written to a file:"
      echo "$filename"_mysql.tar
    else
      url=$(grep -i "/$species/" tmp_$$/ensembl_species.txt)
      name=$(grep -i "/$species/" tmp_$$/ensembl_species.txt | awk -F "/" '{print $(NF-2)}' )
      filename=$(grep -i "/$species/" tmp_$$/ensembl_species.txt | awk -F "/" '{print $(NF)}' )
      echo
      echo "Downloading the genomic sequece of $name"
      echo

      # resolve exact url, because possible http proxy doesn't support wildcards
      wildcard_name=$(basename $url)
      path=$(dirname $url)"/"
      exact_url="$path"$(curl --proxy "" --silent --list-only $path/ | grep -E "$wildcard_name")

      wget -o log "$exact_url" >> /dev/null
      gzipfile=$(ls $filename | grep -i $species)
      echo "Unzipping $gzipfile"
      if [[ $outputmode == "single" ]]
      then
      gunzip $gzipfile
        echo "The results have been written to a file:"
        echo $gzipfile | sed s/".gz"/""/g
      else
        gunzip < $gzipfile >> $outfile
        rm -f $gzipfile
      fi
    fi
  fi


  if [[ $numhits > 1 ]]
  then
    echo " "
    echo "Several genomes match to the given species name"
    echo "Please select a uniqe species name "
    echo ""
    echo "The list of matching species names:"
    grep -i "/$species/" tmp_$$/ensembl_species.txt 

  fi
done

clean_up

