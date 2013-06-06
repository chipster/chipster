#!/bin/bash

tools_path="/opt/chipster/tools/"
export PATH=${PATH}:/opt/chipster4/comp/modules/admin/shell/:/opt/chipster/tools/emboss/bin/

#echo "BWA indexes"
#ls -l $tools_path/bwa_indexes/*.bwt | awk -F "/" '{print "   "$NF}' | sed s/'.fa.bwt'/""/g 
ls -l $tools_path/bwa_indexes/*.bwt | awk -F "/" '{print "   "$NF}' | sed s/'.fa.bwt'/""/g > index_list_$$.bwa.tmp

#echo "Cheking Bowtie indexes"
ls -l $tools_path/bowtie/indexes/*.rev.2.ebwt | awk -F "/" '{print "   "$NF}' | sed s/'.rev.2.ebwt'/""/g >> index_list_$$.bowtie.tmp

#echo "Chreking Bowtie2 indexes"
ls -l $tools_path/bowtie2/indexes/*.rev.1.bt2 | awk -F "/" '{print "   "$NF}' | sed s/'.rev.1.bt2'/""/g >> index_list_$$.bowtie2.tmp

#cheking gtf
ls -l $tools_path/genomes/gtf/*.gtf | grep -v DEXSeq | awk -F "/" '{print "   "$NF}' | sed s/'.gtf'/""/g >> gtf_list_$$.gtf.tmp

#cheking fa
ls -l $tools_path/genomes/fasta/*.fa | awk -F "/" '{print "   "$NF}' | sed s/'\.fa$'/""/g >> index_list_$$.fa.tmp

#cheking fa
ls -l $tools_path/genomes/fasta/nochr/*.fa | awk -F "/" '{print "   "$NF}' | sed s/'\.fa$'/""/g >> index_list_$$.fa_nochr.tmp



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

rm -f index_list_$$.bwa.tmp index_list_$$.bowtie.tmp index_list_$$.bowtie2.tmp index_list_$$.gtf.tmp index_list_$$.fa.tmp index_list_$$.fa_nochr.tmp
