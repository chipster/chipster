#!/bin/bash -f
# A script to write a report about a sequence file or database
# 21.1. 2011 KM

#Syntax:
# infoseq_summary.sh emboss_path sequence_file  
#

##
# sfcheck aliohjelma
##
sfcheck () 
{

if [[ -e $1/infoseq ]]
then
  a="OK"
else
   echo "EMBOSS commands were not found"
  echo " "
  echo "Please provide the path to the EMBOSS executables as the first argument "
  echo "  "
  exit 1
fi


if [[ -e $2 ]] 
then
  a="OK"
else
  echo "Input file $2 not found"
  echo " "
  echo "Please provide a name of a sequence file to be checked as the second argument "
  echo "  "
  exit 1
fi

export LC_ALL=C
export PATH=$1:${PATH}


fformat=$(file $2 |awk '{print $2}')

#bam tiedostyt pitää lukea kokonaan
if [[ $fformat == "gzip" ]]
then
  format=$(infoseq -nowarning -noerror -nodie -only -usa -nohead -filter $2 | awk -F "::" '{print $1}' | uniq | head -1)
else
  format=$(head -800 $2 | infoseq -nowarning -noerror -nodie -only -usa -nohead -filter | awk -F "::" '{print $1}' | uniq | head -1)
  if [[ $format == "" ]]
  then
      format=$(head -99999 $2 | infoseq -nowarning -noerror -nodie -only -usa -nohead -filter | awk -F "::" '{print $1}' | uniq | head -1)
  fi
fi


if [[ $format == "fastq" ]]
then
   q_count=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | wc -c)
   Sa_Il18=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '!"#$%&()*+,-./0123456789:' | wc -c )
   So_1l13_Il15=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d 'KLMNOPQRSTUVWXYZ[\]^_abcdefgh' |wc -c )

   #echo Sanger or Illumina 1.8
   if [[ $q_count != $Sa_Il18 ]]
   then
     Il_18=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '!#$%&()*+,-./0123456789:;<=>?@ABCDEFGHI' | wc -c )
     if [[ $Il_18 -gt 0 ]]
     then
        #echo "Fastq format is: Illumina 1.8+"
        format=("fastq Phred+33 illumina18+")
     else
       #echo "Fastq format is: Sanger or Illumina 1.8+"
       format=("fastq Phred+33")
     fi
  fi

  #echo Solexa or Illumina 13+ or Illumina 15+ 
  if [[ $q_count != $So_1l13_Il15 ]]
  then
    So=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d ';<=>?' | wc -c ) 
    if [[ $q_count -ne $So ]] 
    then
      #echo "Fastq format is: Solexa" 
      format=("fastq Solexa+64")
    else
      So_Il13=$(head -800 $2 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '@A' | wc -c )
      if [[ $q_count -ne $So_Il13 ]] 
         then
         #echo "Fastq format is: Solexa or Illumina 13+" 
         format=("fastq Phred+64")
      else
         #echo "Fastq format is: Solexa or Illumina 13+ or Illumina 15+ " 
        format=("fastq Phred+64")
      fi
    fi 
  fi
fi

if [[ "$format" == "" ]] 
then
  echo "Not an EMBOSS compatible sequence file"
else
  echo $format
fi

}



export LC_ALL=C
export PATH=$1:${PATH}


if [[ "$2" == "" ]]
 then
   echo "Syntax:"
   echo "infoseq_summary emboss_path sequence_file"
   exit 1
fi

if [[ -e $2 ]]
then
   infile=$2
#   printf "%s\t%s\n" " File name:                           " $2
   size=$(du -sh "$2" | awk '{print $1}' )
   printf "%s\t%s\n" " File size:                          " $size
   sformat=$(sfcheck $1 "$2")
   printf "%s\t%s\n" " File format:                         " "$sformat"
   if [[ "$sformat" == "Not an EMBOSS compatible sequence file" ]]
   then
     exit 0
   fi
else 
    echo ' USA definition:                      '"\t""$2"
fi

infoseq -nowarning -nocolumn -delimiter "::"  -nohead "$2" -only -usa -name -type -length -filter | awk -F "::" 'BEGIN { s = 10000000000} { a = a +$NF} \
{ if ( $NF > l) l = $NF } { if ( $NF == l) ln = $3 }  { if ( $NF < s) s = $NF} { if ( $NF == s) sn = $3} {ka = a / NR} \
END { if ( $4 == "N")  print " Sequence type:                       \tNucleotide"} \
END { if ( $4 == "P")  print " Sequence type:                       \tProtein"} \
END { print " Number of sequences:                  \t" NR } \
END { print " Longest (or one of equally long ones): \t" ln "\t" l  } \
END { print " Shortest (or one of equally short ones):\t" sn "\t"s } \
END { print " Average length:                      \t" ka } \
END {print  " Total amount of nucleotides/residues:\t" a } \
END { if ( NF > 5) print " Note: Sequence names seem to contain douple-douple point (::) characters!"}'

exit 0
