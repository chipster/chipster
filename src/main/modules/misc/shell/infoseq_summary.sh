#!/bin/bash -f
# A script to write a report about a sequence file or database
# 21.1. 2011 KM

#Syntax:
# infoseq_summary.sh emboss_path sequence_file  
#

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
   sformat=$(sfcheck.bash $1 "$2")
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
END { print " Number of sequeces:                  \t" NR } \
END { print " Longest (or one of equally long ones): \t" ln "\t" l  } \
END { print " Shortest (or one of equally short ones):\t" sn "\t"s } \
END { print " Average length:                      \t" ka } \
END {print  " Total amount of nucleotides/residues:\t" a } \
END { if ( NF > 5) print " Note: Sequence names seem to contain douple-douple point (::) characters!"}'

exit 0
