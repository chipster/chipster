#!/bin/bash

chipster_path="0"

while [[ $# -ge 1 ]]
do
  case "$1" in
             '-chipster_path')
	      chipster_path="$2"
                shift
                shift
              ;;
              #
  esac
done


if [[ $chipster_path == "0" ]]
then
  echo ""  
  echo "Please define the location of the tools directory of Chipster with option:"
  echo "  -chipster_path /path/to/tools"
  echo
  exit 1
fi 

tools_path="$chipster_path""/tools"
comp_path="$chipster_path""/comp"

export PATH=${PATH}:$chipster_path/comp/modules/admin/shell/:$chipster_path/tools/emboss/bin/

#echo "BWA indexes"
#ls -l $tools_path/bwa_indexes/*.bwt | awk -F "/" '{print "   "$NF}' | sed s/'.fa.bwt'/""/g 
ls -l $tools_path/bwa_indexes/*.bwt | awk -F "/" '{print "   "$NF}' | sed s/'.fa.bwt'/""/g | tr -d " " > index_list_$$.bwa.tmp
grep PARAMETER $comp_path/modules/ngs/R-2.12/bwa* | grep genome: | awk -F "[" '{print $2}' | awk -F "]" '{print $1}' | tr "," "\n" | awk -F : '{print $1}' | sort | uniq | tr -d  " " |  sed s/'\.fa$'/""/g > genomes_in_bwa_interface.$$.tmp

echo "--------------------------------------------------------------------------------------------------------------------------"
echo "BWA index summary"
cat index_list_$$.bwa.tmp genomes_in_bwa_interface.$$.tmp genomes_in_bwa_interface.$$.tmp | sort | uniq -c | \
awk '{ if ($1 == "1" ) printf "%60s\tIndexed but not in use\n", $2 }{if ($1 == "2" ) printf "%60s\tIndex missing!\n", $2}{if ( $1 == "3" ) printf "%60s\tOK\n", $2}'
echo "---------------------------------------------------------------------------------------------------------------------------"


#echo "Cheking Bowtie indexes"
ls -l $tools_path/bowtie/indexes/*.rev.2.ebwt | awk -F "/" '{print "   "$NF}' | sed s/'.rev.2.ebwt'/""/g | tr -d " " >> index_list_$$.bowtie.tmp

grep PARAMETER $comp_path/modules/ngs/R-2.12/bowtie[.-]* | grep genome: | awk -F "[" '{print $2}' | awk -F "]" '{print $1}' | tr "," "\n" | awk -F : '{print $1}' | sort | uniq | tr -d  " " |  sed s/'\.fa$'/""/g > genomes_in_bowtie_interface.$$.tmp
#echo "--------------------------------------------------------------------------------------------------------------------------"
echo "Bowtie index summary"
cat index_list_$$.bowtie.tmp genomes_in_bowtie_interface.$$.tmp genomes_in_bowtie_interface.$$.tmp | sort | uniq -c | \
awk '{ if ($1 == "1" ) printf "%60s\tIndexed but not in use\n", $2 }{if ($1 == "2" ) printf "%60s\tIndex missing!\n", $2}{if ( $1 == "3" ) printf "%60s\tOK\n", $2}'
echo "---------------------------------------------------------------------------------------------------------------------------"

#echo "Chreking Bowtie2 indexes"
ls -l $tools_path/bowtie2/indexes/*.rev.1.bt2 | awk -F "/" '{print "   "$NF}' | sed s/'.rev.1.bt2'/""/g | tr -d " " >> index_list_$$.bowtie2.tmp

grep PARAMETER $comp_path/modules/ngs/R-2.12/bowtie2* | grep genome: | awk -F "[" '{print $2}' | awk -F "]" '{print $1}' | tr "," "\n" | awk -F : '{print $1}' | sort | uniq | tr -d  " " |  sed s/'\.fa$'/""/g > genomes_in_bowtie2_interface.$$.tmp
#echo "--------------------------------------------------------------------------------------------------------------------------"
echo "Bowtie2 index summary"
cat index_list_$$.bowtie2.tmp genomes_in_bowtie2_interface.$$.tmp genomes_in_bowtie2_interface.$$.tmp | sort | uniq -c | \
awk '{ if ($1 == "1" ) printf "%60s\tIndexed but not in use\n", $2 }{if ($1 == "2" ) printf "%60s\tIndex missing!\n", $2}{if ( $1 == "3" ) printf "%60s\tOK\n", $2}'
echo "---------------------------------------------------------------------------------------------------------------------------"

#cheking gtf
ls -l $tools_path/genomes/gtf/*.gtf | grep -v DEXSeq | awk -F "/" '{print "   "$NF}' | sed s/'.gtf'/""/g | tr -d " " >> gtf_list_$$.gtf.tmp

#cheking fa
ls -l $tools_path/genomes/fasta/*.fa | awk -F "/" '{print "   "$NF}' | sed s/'\.fa$'/""/g | tr -d " " >> index_list_$$.fa.tmp

#cheking fa
ls -l $tools_path/genomes/fasta/nochr/*.fa | awk -F "/" '{print "   "$NF}' | sed s/'\.fa$'/""/g | tr -d " " >> index_list_$$.fa_nochr.tmp


echo "Summary of genome indexes and files"
echo
printf "%60s %7s %7s %7s %7s %7s %7s \n" Species BWA Bowtie Bowtie2 gtf fasta fasta_nochr
echo "----------------------------------------------------------------------------------------------------------"


for ind in $(sort index_list_$$.*.tmp | uniq)
do
 bwa_m=$(grep -c $ind index_list_$$.bwa.tmp)
 if [[ bwa_m -eq 1 ]]
 then
    bwa_m="OK"
 else
    bwa_m="-"
 fi

 bt_m=$(grep -c $ind index_list_$$.bowtie.tmp)
 if [[ bt_m -eq 1 ]]
 then
    bt_m="OK"
 else
    bt_m="-"
 fi

 bt2_m=$(grep -c $ind index_list_$$.bowtie2.tmp)
 if [[ bt2_m -eq 1 ]]
 then
    bt2_m="OK"
 else
    bt2_m="-"
 fi

 ind2=$(echo $ind | sed s/"\.dna\.toplevel"/""/g )  
 gtf_m=$(grep -c $ind2 gtf_list_$$.gtf.tmp)
 if [[ gtf_m -eq 1 ]]
 then
    gtf_m="OK"
 else
    gtf_m="-"
 fi

 fa_m=$(grep -c $ind index_list_$$.fa.tmp)
 if [[ fa_m -eq 1 ]]
 then
    fa_m="OK"
 else
    fa_m="-"
 fi

 fa_nochr_m=$(grep -c $ind index_list_$$.fa_nochr.tmp)
 if [[ fa_nochr_m -eq 1 ]]
 then
    fa_nochr_m="OK"
 else
    fa_nochr_m="-"
 fi

 printf "%60s %7s %7s %7s %7s %7s %7s \n" $ind $bwa_m $bt_m $bt2_m $gtf_m $fa_m $fa_nochr_m

# echo $ind $bwa_m $bt_m $bt2_m
done

rm -f *_$$.*.tmp 



