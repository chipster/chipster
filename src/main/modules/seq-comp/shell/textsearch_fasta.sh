#!/bin/bash 
# A script to combine textsearch and seqret
# 29.10. 2013 KM

#setenv LC_ALL C
#setenv PLPLOT_LIB /v/linux26_x86_64/appl/molbio/emboss/6.3.1/share/EMBOSS
#setenv PATH /v/linux26_x86_64/appl/molbio/emboss/6.3.1/bin:$PATH
#setenv LD_LIBRARY_PATH /v/linux26_x86_64/appl/molbio/emboss/libharu/2.0.8/lib
emboss_path=("/appl/bio/emboss/6.5.7_gnu/bin")
osformat=("fasta")
casesensitive=("N")

if [[ "$1" == "" ]] 
then
   echo "Syntax:  textsearch_fasta -emboss_path path -sequence seqfile.fasta -pattern \"pattern1|pattern2\" -outfile results.fasta"
   exit
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
              '-pattern')
                pattern=($2)
               shift
               shift
               ;;
              '-casesensitive')                    # Data tyoe
                casesensitive=($2)                        
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
               shift  
               exit
               ;;

    esac
done
export PATH=${emboss_path}:${PATH}

textsearch -sequence $sequence -pattern "$pattern" -casesensitive $casesensitive -outfile ${outfile}_"$$"_list -only -usa -auto 2> /dev/null
test=$(head  ${outfile}_"$$"_list | wc -l)

if [[ $test -eq 0 ]] 
then
   echo "No sequences matching search pattern: $pattern were found" > $outfile
   rm -rf "$outfile"_"$$"_list
   exit 
else
   seqret @"$outfile"_"$$"_list -outseq $outfile -auto
fi

rm -rf "$outfile"_"$$"_list

exit 1
