#!/bin/bash
# A script to filter protein sequece sets
# 23.1. 2014 KM

export LC_ALL=C
#setenv PLPLOT_LIB /v/linux26_x86_64/appl/molbio/emboss/6.5.7/share/EMBOSS
#export PATH=/appl/molbio/emboss/6.5.7_gnu/bin:${PATH}
#setenv LD_LIBRARY_PATH /v/linux26_x86_64/appl/molbio/emboss/libharu/2.0.8/lib
osformat=("fasta")
casesensitive=("N")
outfile=("stdout")
text=(0)
fuzzpro=(0)
minlength=(-1)
maxlength=(-1)
minmass=(-1)
maxmass=(-1)
emboss_path=(x)

if [[ "$1" == "" ]]
then

   echo "protfilter command line options:"
   echo " -sequence, -desctext, -casesensitive, -pattern, -maxlength, -minlength, -minmass, -maxmass, -outfile  "
   exit 1
fi
while [[ $# -ge 1 ]]
do
  case "$1" in
              '-emboss_path')
               emboss_path=($2)
               shift
               shift
               ;;  
               '-sequence')
               sequence=($2)
               shift
               shift
               ;;
               '-desctext')
                desctext=($2)
                text=(1)
                shift
                shift
              ;;
              '-casesensitive')                   
                casesensitive=($2)                        
                shift
                shift
              ;;
              '-pattern')
                fuzzpro=(1)
                pattern=("$2")
                shift
                shift
              ;;
              '-maxlength')
                length=(1)
                maxlength=($2)
                shift
                shift
              ;;
              '-minlength')
                length=(1)
                minlength=($2)
                shift
                shift
              ;;
              '-minmass')
                mass=(1)
                minmass=($2)
                shift 
                shift 
              ;;
              '-maxmass')
                mass=(1)
                maxmass=($2)
                shift
                shift
              ;;
              '-outfile')
                outfile=($2) 
                  if [[ -e "$outfile" ]] 
                  then
                     echo "Result file: $outfile already exists!"
                     exit 1
                  fi
                  shift
                  shift 
              ;;
              *)
                echo "Unknown argument or parameter: $1"
                shift argv 
                exit 1
              ;;

    esac
done

if [[ $emboss_path == x ]] 
then 
   emboss_path=$(which seqret | sed s/"\/seqret"/""/g)
   if [[ $emboss_path == "" ]]
   then 
      echo emboss commands not found
      exit 1
   fi
fi

$emboss_path/infoseq -nohead -only -usa $sequence -outfile infile_"$$"_list -filter 
nseq=$(wc -l infile_"$$"_list | awk '{print $1}')
echo "The input sequece set $sequence contains $nseq sequcenses"


if [[ $fuzzpro -eq 1 ]] 
then
  $emboss_path/fuzzpro -sequence @infile_"$$"_list -pattern "$pattern" -outfile outfile_"$$"_fuzzpro -rformat listfile -auto
  awk -F "[" '{print $1}' outfile_"$$"_fuzzpro | grep "::" | grep -v "^#" | sort | uniq > outfile_"$$"_list
  rm -f outfile_"$$"_fuzzpro
  nseq=$(wc -l outfile_"$$"_list |awk '{print $1}') 
  if [[ $nseq -eq 0 ]] 
  then
     echo "No sequences matching search pattern: $pattern were found"
     echo "No sequences matching search pattern: $pattern were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
     echo "Sequences matching search pattern: $pattern : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list
fi

if [[ $maxlength -gt -1 ]]
then 
   $emboss_path/infoseq @infile_"$$"_list -nohead -only -usa -length -filter | awk '{if ( $2 <='"$maxlength"' ) print $1}' > outfile_"$$"_list
   nseq=$(wc -l outfile_"$$"_list |awk '{print $1}') 
   if [[ $nseq -eq 0 ]] 
   then
     echo "No sequences matching condition maxlength=$maxlength were found"
     echo "No sequences matching condition maxlength=$maxlength were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
     echo "Sequences matching condition maxlength=$maxlength : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list
fi


if [[ $minlength -gt -1 ]] 
then
   $emboss_path/infoseq @infile_"$$"_list -nohead -only -usa -length -filter | awk '{if ( $2 >='"$minlength"' ) print $1}' > outfile_"$$"_list
  nseq=$(wc -l outfile_"$$"_list |awk '{print $1}')
  if [[ $nseq -eq 0 ]] 
  then
     echo "No sequences matching condition minlengh=$minlength were found"
     echo "No sequences matching condition minlengh=$minlength were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
     echo "Sequences matching condition minlengh=$minlength : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list
fi

if [[ $minmass -gt -1 ]] 
then
  for seq in $(cat infile_"$$"_list)
  do
     mass=$($emboss_path/pepstats $seq -filter | grep "Molecular weight" | awk '{ print $4 }'| sed s/"\."/""/g )
     (( minmass100 = minmass * 100 ))
     if [[ $mass -ge $minmass100 ]] 
     then
         echo $seq >>! outfile_"$$"_list
     fi
  done
  if [[ ! -e outfile_"$$"_list  ]] 
  then
     echo "No sequences matching condition minmass= $minmass were found"
     echo "No sequences matching condition minmass= $minmass were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
      nseq=$(wc -l outfile_"$$"_list |awk '{print $1}' )
      echo "Sequences matching condition minmass= $minmass : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list
fi

if [[ $maxmass -gt -1 ]] 
then
  for seq in $(cat infile_"$$"_list)
  do
     mass=$($emboss_path/pepstats $seq -filter | grep "Molecular weight" | awk '{ print $4 }'| sed s/"\."/""/g )
     (( maxmass100 = maxmass * 100 ))
     if [[ $mass -le $maxmass100 ]] 
     then
         echo $seq >>! outfile_"$$"_list
     fi
  done
  if [[ ! -e outfile_"$$"_list  ]] 
  then
     echo "No sequences matching condition maxmass= $maxmass were found"
     echo "No sequences matching condition maxmass= $maxmass were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
      nseq=$(wc -l outfile_"$$"_list |awk '{print $1}' )
      echo "Sequences matching condition maxmass= $minmass : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list
fi  


if [[ $text -eq 1 ]] 
then
  $emboss_path/textsearch  -sequence @infile_"$$"_list -pattern "$desctext" -casesensitive $casesensitive -outfile outfile_"$$"_list -only -usa -auto
  nseq=$(wc -l outfile_"$$"_list |awk '{print $1}') 
  if [[ $nseq == 0 ]] 
  then
     echo "No sequence descriptions matching search pattern: $desctext were found" 
     echo "No sequence descriptions matching search pattern: $desctext were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  else
     echo "Sequence descriptions matching search pattern: $desctext : $nseq"
  fi
  rm -f infile_"$$"_list
  mv outfile_"$$"_list infile_"$$"_list 
fi


 $emboss_path/seqret @infile_"$$"_list -outseq $outfile -auto
rm -rf infile_"$$"_list
