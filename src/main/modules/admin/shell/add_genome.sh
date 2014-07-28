#!/bin/bash 

#export PATH=${PATH}:/opt/chipster4/comp/modules/admin/shell/:/opt/chipster/tools/emboss/bin/:/opt/chipster/tools/samtools/:/opt/chipster/tools/tabix

add_genome_help ()
{
cat <<EOF
------------------------------------------------------------------------------------
add_genome.sh script adds new genome to the Chipster.
Syntax:

  add_genome.sh -chipster_path chipster_tools_path species_name

The command abowe retrieves form Ensembl or Ensembl genomes database the sequence 
(topleve.fasta) and description (gtf) files for the given species. After that
BWA, Bowtie and Bowtie2 indexes are caluculated for the genomic sequence. 
The MySQL files are retrieved too.

The gtf and MySQL files are processed to create the files used by the Genome Browser.

Other options:

  -only_bwa       Calculate only bwa indexes
  -only_bowtie    Calculate only bowtie indexes
  -only_bowtie2   Calculate only bowtie2 indexes

------------------------------------------------------------------------------------


EOF
}

####
#Process gtf aliohjelma
####
process_gtf () # parameters 1:url
{
	FILE_BODY=$(basename $1 .gtf)

	#generate list of chromosomes of genes
	#Read file  		Take only chr and name columns     	Filter out other names    	Remove duplicates   Replace useless chars with tab Or remove        And write to file
	cat "$FILE_BODY.gtf" | 	cut -f 1,9 --output-delimiter=';' | 	cut -d ';' -f 1,5      | 	uniq              | sed -e 's/; gene_name "/	/' | sed -e 's/\"//' > "$FILE_BODY.gene.tsv"


	#tabix installation folder hast to be in $PATH to find bgzip and tabix programs
	#don't exit even if grep exits with error (when there isn't any comments)
	set +e
	grep "^#" "$FILE_BODY.gtf" > "$FILE_BODY-1.gtf"
	set -e

	grep -v "^#" "$FILE_BODY.gtf" >> "$FILE_BODY-1.gtf"
	cat $FILE_BODY-1.gtf | sort -k1,1 -k4,4n -t "	" > "$FILE_BODY-sorted.gtf"		
	cat "$FILE_BODY-sorted.gtf" | bgzip > "$FILE_BODY.tabix.gtf.gz"

	#rm "$FILE_BODY.gtf"
	rm "$FILE_BODY-1.gtf"
	rm "$FILE_BODY-sorted.gtf"

	#generate index
	tabix -p gff "$FILE_BODY.tabix.gtf.gz";

#Result files	
#"$FILE_BODY.gene.tsv"
#"$FILE_BODY.tabix.gtf.gz"
#"$FILE_BODY.tabix.gtf.gz.tbi"	
}


# process ensembl mysql files
ensembl_mysql () # parameters 1:url 2:new-name
{

	echo Prepare chromosome name and identifier mapping
        set +e

	# search for chomosome (or group) coordinate systems (group for stickleback)
        grep "chromosome\|group" $1/coord_system.txt > coord_system-chr.txt

	# join requires sorted input
	LANG=en_EN sort -k 1 coord_system-chr.txt > coord_system-sorted.txt
	LANG=en_EN sort -k 3 "$1/seq_region.txt" > seq_region-sorted.txt

	# join chromosome names and seq_region identifiers to create map of chromosome identifiers
	LANG=en_EN join -t '	' -1 1 -2 3 coord_system-sorted.txt seq_region-sorted.txt > chr_map-join.txt
	
	# remove extra columns
	cat chr_map-join.txt | cut -f 7,8 > chr_map.txt

	# join requires sorted input
	LANG=en_EN sort -k 1 chr_map.txt > chr_map-sorted.txt
	
	# clean
	rm "$1/coord_system.txt" "$1/seq_region.txt"
	rm coord_system-chr.txt coord_system-sorted.txt seq_region-sorted.txt chr_map-join.txt chr_map.txt

	
	# Low complexity region data

	# remove extra columns
	#cat repeat-masker-join.txt | cut -f 3,4,5 > repeat-masker.txt		
	cat "$1/repeat_feature.txt" | cut -d '	' -f 2,3,4 > repeat-masker.txt
		
	# join requires sorted input
	LANG=en_EN sort -k 1 repeat-masker.txt > repeat-masker-sorted.txt

	# join chromosome identifiers and the data
	LANG=en_EN join -t '	' -1 1 -2 1 chr_map-sorted.txt repeat-masker-sorted.txt > repeat-join.txt

	# remove extra columns
	# FIXME Ensembl uses 1-based coordinates, whereas standard bed file must have 0-based coordinates
	cat repeat-join.txt | cut -d '	' -f 2,3,4 > repeat.bed
	
	# bed to tabix
	cat repeat.bed | sort -k1,1 -k2,2n > repeat-sorted.bed		
	cat repeat-sorted.bed | bgzip > "$1/repeat-tabix.bed.gz"

	#generate index
	tabix -p bed "$1/repeat-tabix.bed.gz"

	# clean
	#rm "$1/analysis.txt" "$1/repeat_feature.txt"
	#rm repeat-masker-row.txt repeat-masker-id.txt repeat-sorted.txt repeat-masker-join.txt
	rm repeat-masker.txt  repeat-masker-sorted.txt repeat-join.txt
	rm repeat.bed repeat-sorted.bed


	# Cytoband data

	# sort the cytoband data
        cp karyotype.txt cytoband-tmp.txt
	LANG=en_EN sort -k 2 "$1/cytoband-tmp.txt" > cytoband-sorted.txt

	# join chromosome identifiers and the cytobands
	LANG=en_EN join -t '	' -1 1 -2 2 chr_map-sorted.txt cytoband-sorted.txt > cytoband-join.txt

	# remove extra columns
	cat cytoband-join.txt | cut -d '	' -f 2,3,4,5,6,7 > "$1/cytoband-chr.txt"

	# clean
	rm repeat_feature.txt karyotype.txt
	rm cytoband-sorted.txt cytoband-join.txt "$1/cytoband-tmp.txt"
	rm chr_map-sorted.txt



	#result files: "$2cytoband-chr.txt", "$2repeat-tabix.bed.gz" and "$2repeat-tabix.bed.gz.tbi"
}

chipster_path="0"
ensembl=0
fasta=0
gtf=0
text=0
length=0
karyotype=0
clean=0
version="0.0"
ensembl_version="0"
location=$(pwd)
INDEX_BWA=1
INDEX_BOWTIE=1
INDEX_BOWTIE2=1
INDEX_TOPHAT2=1

while [[ $# -ge 1 ]]
do
  case "$1" in
              '-chipster_path')
	      chipster_path="$2"
                shift
                shift
              ;;
              '-genomes_path')
              genomes_path="$2"
                shift
                shift
              ;;
              #
              '-species')
                species="$2"
                ensembl=0
                shift 
                shift 
              ;;
              '-fasta')
                genome_fasta="$2"
                fasta=1
                shift 
                shift 
              ;;
              '-gtf')
                genome_gtf="$2"
                gtf=1
                shift 
                shift 
              ;;
              '-version')
                ensembl_version="$2"
                shift 
                shift 
              ;;
              '-descfilter')
                desctext="$2"
                text="1"
                shift 
                shift 
              ;;
              '-minlength')
                length=1
                minlength="$2"
                shift
                shift
              ;;
              '-karyotype')
                karyotype=1
                shift
              ;;
              '-clean')
              clean=1
                shift
              ;;
              '-only_bwa')
                INDEX_BOWTIE=0
                INDEX_BOWTIE2=0
                shift
              ;;
	      '-only_bowtie')
                INDEX_BOWTIE2=0
                INDEX_BWA=0
                shift
              ;;
	      '-only_bowtie2')
                INDEX_BOWTIE=0
                INDEX_BWA=0
                shift
              ;;
	      '-index')
                INDEX_BOWTIE=1
                INDEX_BOWTIE2=1
                INDEX_BWA=1
                shift
              ;;
	      '-no_index')
                INDEX_BOWTIE=0
                INDEX_BOWTIE2=0
                INDEX_BWA=0
                shift
              ;;
              '-help')
              add_genome_help
              exit 0
              ;;
              *)
                species="$1"
                ensembl=1
                shift 
              ;; 

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

species=$(echo $species | sed s/" "/"_"/g )

tools_path=${chipster_path}/tools
comp_path=${chipster_path}/comp
if [[ -z "$genomes_path" ]]
then
  genomes_path=${tools_path}/genomes
fi
index_path=${genomes_path}/indexes
tmp_path=${genomes_path}/tmp/${species}_$$ # process id

if [[ $clean -eq 1 ]]
then
  echo "debug: removing $genomes_path/genomebrowser/$species"
  echo "debug: $(ls -lah $genomes_path/genomebrowser/$species)"
  rm -rf ${genomes_path}/tmp/${species}_* $genomes_path/genomebrowser/$species $genomes_path/fasta/$species* $genomes_path/gtf/$species* $index_path/bowtie/$species* $index_path/bowtie2/$species* $index_path/bwa/$species* $index_path/tophat2/$species*
  echo "debug: $(ls -lah $genomes_path/genomebrowser/$species)"
fi

# these will fail occasionally when this script is run in parallel
if [[ ! -e ${genomes_path} ]]
then
  mkdir --parents ${genomes_path}
fi

if [[ ! -e ${index_path} ]]
then
  mkdir --parents ${index_path}
fi

if [[ -e ${tmp_path} ]]
then
  echo "Work directory ${tmp_path} exists already. Please terminate all scripts using this directory and delete it."
  exit 1
else
  mkdir --parents ${tmp_path}
fi

# exit if there is an error
set -e

export PATH=${PATH}:$comp_path/modules/admin/shell/:$tools_path/emboss/bin/:$tools_path/samtools/:$tools_path/tabix/:$tools_path/bowtie2/
echo $PATH

echo "Installing genomes to: $genomes_path"

##
#  Retrieve the fasta file
##

#reading the data from ensembl

echo $ensembl

cd $tmp_path

if [[ $ensembl -eq 1 ]]
then
  echo "Retrieving and indexing genome sequence for $species"

  if [[ $ensembl_version -eq 0 ]]
  then
    #echo ensemblfetch.sh $species
    genome_fasta=$(ensemblfetch.sh $species | tail -1)
  else
    #echo ensemblfetch.sh $species -version $ensembl_version
    genome_fasta=$(ensemblfetch.sh $species -version $ensembl_version | tail -1)
  fi

  genome_name=$(basename $genome_fasta .dna.toplevel.fa)

  if [[ $genome_fasta == "--------------------------------------------------------------------------------" ]]
  then
    echo Species $species was not found from the Ensembl database.
    exit 1
  fi

  #If karyotype option is used, pick only chromosomes
  if [[ $karyotype -eq 1 ]]
  then
    echo "Filtering karyotype chromosomes..."
    curl -s "http://beta.rest.ensembl.org/assembly/info/$species?content-type=text/xml" | grep "<karyotype>" | sed "s/<karyotype>//g" | sed "s/<\/karyotype>//g" | sed "s/ *//g" > karyotype.tmp
    curl -s "http://beta.rest.ensemblgenomes.org/assembly/info/$species?content-type=text/xml" | grep "<karyotype>" | sed "s/<karyotype>//g" | sed "s/<\/karyotype>//g" | sed "s/ *//g" >> karyotype.tmp
    num_chrs=$(cat karyotype.tmp | wc -l )
    if [[ $num_chrs > 0 ]]
    then
       mv $genome_fasta toplevel.fasta
       awk '{print "toplevel.fasta:"$1}' karyotype.tmp  >  karyotype.list
       seqret @karyotype.list $genome_name.fa
       rm -f toplevel.fasta karyotype.list
    else
       echo "no karyotype chromosomes, keeping all"
       mv $genome_fasta $genome_name.fa
    fi
    genome_fasta=$genome_name.fa

    rm -f karyotype.tmp
  fi

  echo "$genome_fasta Downloaded"

  echo "Creating fasta index..."

  #Calculate index file for the fasta file
  ${tools_path}/samtools/samtools faidx ${genome_fasta}

  echo "Content:"
  infoseq $genome_fasta -auto

  #version=$(echo $genome_fasta | awk -F "." '{print $2}')
  version=$(echo $genome_fasta | awk -F "." '{for (i=2; i<=NF; i++) printf $i"." }' | awk -F ".fa" '{print $1}')
fi


if [[ $fasta -eq 1 ]]
then
  #With this option, do not fetch fasta from ensembl but use given file instead
  cp $genome_fasta ${tmp_path}/
fi


if [[ "$text" == "1" ]]
  then
  textsearch  -sequence $genome_fasta -pattern "$desctext" -outfile outfile_"$$"_list -only -usa -auto >> /dev/null
  test=$(head outfile_"$$"_list | wc -l) 
  if [[ $test -eq 0 ]] 
  then 
     echo "No sequences matching search pattern: $desctext were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  fi
  seqret @outfile_"$$"_list -outseq ${genome_fasta}.filtered -auto
  rm -f $genome_fasta
  mv ${genome_fasta}.filtered $genome_fasta
fi

###
#  get the gtf file
###


if [[ $ensembl -eq 1 ]]
then

  #which ensemblfetch.sh

  if [[ $ensembl_version -eq 0 ]]
  then
    #echo ensemblfetch.sh -type gtf $species
    genome_gtf=$(ensemblfetch.sh -type gtf $species | tail -1 )
  else
    #echo ensemblfetch.sh -type gtf $species -version $ensembl_version
    genome_gtf=$(ensemblfetch.sh -type gtf $species -version $ensembl_version | tail -1 )
  fi
  
  genome_gtf_name=$(basename $genome_gtf .gtf)
  echo "executing: python ${tools_path}/dexseq-exoncounts/dexseq_prepare_annotation.py $genome_gtf $genome_name.DEXSeq.gtf "
  python ${tools_path}/dexseq-exoncounts/dexseq_prepare_annotation.py $genome_gtf $genome_name.DEXSeq.gtf 

  process_gtf $genome_gtf 

fi

if [[ $gtf -eq 1 ]]
then  
  cp $location/$genome_gtf ${genomes_path}/gtf/
  python ${tools_path}/dexseq-exoncounts/dexseq_prepare_annotation.py $genome_gtf $genome_name.DEXSeq.gtf 
  process_gtf $genome_gtf
fi


###
#  get the mysql files
###
echo "downloading mysql files"
if [[ $ensembl -eq 1 ]]
then

  if [[ $ensembl_version -eq 0 ]]
  then
    mysql_files=$(ensemblfetch.sh -type mysql $species | tail -1)
  else
    mysql_files=$(ensemblfetch.sh -type mysql $species -version $ensembl_version | tail -1 )
  fi

  testtsring=$(echo $mysql_files | wc -c)
  if [[ $testtsring -gt 50 ]]
  then
    echo "No MySQL data for $species was found from the Ensembl database."
  else    
    mysql_dir=$(basename $mysql_files .tar)
    tar xf $mysql_files
    gunzip $mysql_dir/*.gz
    mv $mysql_dir/* .
    rm -f $mysql_files
    rmdir $mysql_dir
    ensembl_mysql .
  fi
fi


##genome.yaml

echo "species: $species" > genome.yaml
echo "version: ($version)" >> genome.yaml
echo "ensemblBrowserUrl: " >> genome.yaml
echo "ucscBrowserUrl: " >> genome.yaml
echo "sortId: other" >> genome.yaml


###
#  Move files to right places
###

gb_path=$genomes_path/genomebrowser/$species/$version

if [[ ! -e ${genomes_path}/fasta ]]
then
  mkdir --parents ${genomes_path}/fasta
fi

if [[ ! -e ${genomes_path}/gtf ]]
then
   mkdir --parents ${genomes_path}/gtf
fi

if [[ ! -e $gb_path ]]
then
  mkdir --parents $gb_path
fi

mv ${genome_fasta} ${genomes_path}/fasta/
mv ${genome_fasta}.fai ${genomes_path}/fasta/
mv ${genome_gtf} ${genomes_path}/gtf/
mv ${genome_name}.DEXSeq.gtf ${genomes_path}/gtf/

ln -s ../../../fasta/$genome_fasta 	$gb_path/$genome_fasta
ln -s ../../../fasta/$genome_fasta.fai 	$gb_path/$genome_fasta.fai

# we don't have mysql data for all genomes
set +e
mv cytoband-chr.txt			$gb_path/
mv genome.yaml				$gb_path/
mv repeat-tabix.bed.gz			$gb_path/
mv repeat-tabix.bed.gz.tbi		$gb_path/
mv ${genome_name}.gene.tsv 		$gb_path/
set -e
mv ${genome_name}.tabix.gtf.gz 		$gb_path/
mv ${genome_name}.tabix.gtf.gz.tbi	$gb_path/




##
#Check if the aligner indexes have been already created
##

#size=$(ls -l $genomes_path/fasta/$genome_fasta | awk '{print $5}')
#checksum=$(md5sum $genomes_path/fasta/$genome_fasta | awk '{print $1}')

#look for matching size and md5sum

#genome_check=$(grep -h $size $tools_path/genomes/genome_list | grep $checksum | awk '{print $1}' | tail -1)


#if [ ! $genome_check == "" ]; then
#  echo "File $genome_fasta has alredy been indexed"
#  exit 0
#fi

#make bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then

  if [[ ! -e ${index_path}/bwa ]]
  then
    mkdir ${index_path}/bwa
  fi

  echo "-------------------------------------------------------"
  echo "Calculating BWA indexes for $genome_fasta"
  cd $index_path/bwa
  ln -s ../../fasta/$genome_fasta $genome_name.fa

  { time $tools_path/bwa/bwa index -p $genome_name  $genome_name.fa ; } &> $tmp_path/bwa.log &

else
  echo "Skipping BWA indexing"
fi


#make bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  if [[ ! -e ${index_path}/bowtie ]]
  then
    mkdir ${index_path}/bowtie
  fi

  echo "-------------------------------------------------------"
  echo "Calculating Bowtie indexes for $genome_fasta"
  cd $index_path/bowtie
  ln -s ../../fasta/$genome_fasta $genome_name.fa
  { time $tools_path/bowtie/bowtie-build $genome_name.fa $genome_name ; } &> $tmp_path/bowtie.log &
else
    echo "Skipping Bowtie indexing"
fi

#make bowtie2_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  if [[ ! -e ${index_path}/bowtie2 ]]
  then
    mkdir ${index_path}/bowtie2
  fi

  echo "-------------------------------------------------------"
  echo "Calculating Bowtie2 indexes for $genome_fasta"

  cd $index_path/bowtie2
  ln -s ../../fasta/$genome_fasta $index_path/bowtie2/$genome_name.fa
  
  { time $tools_path/bowtie2/bowtie2-build $genome_name.fa $genome_name ; } &> $tmp_path/bowtie2.log &

else
    echo "Skipping Bowtie2 indexing"
fi

# wait for all background tasks to complete (bwa, bowtie and bowtie2 indexing)
wait

  #check bowtie indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  cd $index_path/bowtie
  n_genome=$(grep -c "^>" $genome_name.fa)
  n_index=$($tools_path/bowtie/bowtie-inspect -n $genome_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie indexing of genome $genome_name failed"
     exit 1
  else
    echo "Bowtie index of genome $genome_name OK"
  fi
fi

  #check bowtie2 indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  cd $index_path/bowtie2
  n_genome=$(grep -c "^>" $genome_name.fa)
  n_index=$($tools_path/bowtie2/bowtie2-inspect -n $genome_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie2 indexing of genome $genome_name failed"
     exit 1
  else
    echo "Bowtie2 index of genome $genome_name OK"

    if [[ $INDEX_TOPHAT2 -eq 1 ]]
    then
       if [[ ! -e ${index_path}/tophat2 ]]
       then
          mkdir ${index_path}/tophat2
       fi

#       ln -s ../../fasta/$genome_fasta $index_path/tophat2/$genome_name.fa
       
       echo "Building TopHat2 transcriptome index using $genome_gtf"		

       { time $tools_path/tophat2/tophat -G  $genomes_path/gtf/$genome_gtf --transcriptome-index $index_path/tophat2/$genome_name $index_path/bowtie2/$genome_name ; } &> $tmp_path/tophat2.log
	
	set +e
	# no idea why first attempt sometimes fails
	rm -rf $index_path/bowtie2/tophat_out
	rm -rf $index_path/bowtie2/tophat_out
	set -e
    fi
  fi
fi

cd $genomes_path
find genomebrowser/$species/$version/* 	>  $tmp_path/$genome_name.gb.files
find fasta/$genome_name* 		>> $tmp_path/$genome_name.fa.files
find gtf/$genome_name* 			>> $tmp_path/$genome_name.gtf.files

if [[ $INDEX_BOWTIE -eq 1 ]]
then
  find indexes/bowtie/$genome_name* 	>> $tmp_path/$genome_name.bowtie.files
fi

if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  find indexes/bowtie2/$genome_name* 	>> $tmp_path/$genome_name.bowtie2.files
  find indexes/tophat2/$genome_name* 	>> $tmp_path/$genome_name.tophat2.files
fi

if [[ $INDEX_BWA -eq 1 ]]
then
  find indexes/bwa/$genome_name* 		>> $tmp_path/$genome_name.bwa.files
fi

echo ""
echo "---------------------------------------------------------------"

