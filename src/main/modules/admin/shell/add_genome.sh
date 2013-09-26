#!/bin/bash 

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

		# Download database dump files

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


		# sort the actual data
                cp karyotype.txt cytoband-tmp.txt
		LANG=en_EN sort -k 2 "$1/cytoband-tmp.txt" > cytoband-sorted.txt

		# join chromosome identifiers and the data
		LANG=en_EN join -t '	' -1 1 -2 2 chr_map-sorted.txt cytoband-sorted.txt > cytoband-join.txt

		# remove extra columns
		cat cytoband-join.txt | cut -d '	' -f 2,3,4,5,6,7 > "$1/cytoband-chr.txt"

		# clean
		rm cytoband-sorted.txt cytoband-join.txt "$1/cytoband-tmp.txt"
		rm chr_map-sorted.txt

	#result files: "$2cytoband-chr.txt", "$2repeat-tabix.bed.gz" and "$2repeat-tabix.bed.gz.tbi"
}

process_mysql_files ()
{

		# Download database dump files
	
#		download_and_rename "$1seq_region.txt.gz" "$2seq_region.txt"
#		download_and_rename "$1coord_system.txt.gz" "$2coord_system.txt"
#		download_and_rename "$1karyotype.txt.gz" "$2cytoband-tmp.txt" 
#		download_and_rename "$1repeat_feature.txt.gz" "$2repeat_feature.txt"

                cp karyotype.txt cytoband-tmp.txt

		# Prepare chromosome name and identifier mapping

		# search for chomosome (or group) coordinate systems (group for stickleback)
		cat "coord_system.txt" | grep "chromosome\|group" > coord_system-chr.txt

		# join requires sorted input
		LANG=en_EN sort -k 1 coord_system-chr.txt > coord_system-sorted.txt
		LANG=en_EN sort -k 3 seq_region.txt > seq_region-sorted.txt

		# join chromosome names and seq_region identifiers to create map of chromosome identifiers
		LANG=en_EN join -t '	' -1 1 -2 3 coord_system-sorted.txt seq_region-sorted.txt > chr_map-join.txt
	
		# remove extra columns
		cat chr_map-join.txt | cut -f 7,8 > chr_map.txt

		# join requires sorted input
		LANG=en_EN sort -k 1 chr_map.txt > chr_map-sorted.txt
	
		# clean
		rm coord_system.txt seq_region.txt
		rm coord_system-chr.txt coord_system-sorted.txt seq_region-sorted.txt chr_map-join.txt chr_map.txt

	
		# Low complexity region data

		# remove extra columns
		#cat repeat-masker-join.txt | cut -f 3,4,5 > repeat-masker.txt		
		cat repeat_feature.txt | cut -d '	' -f 2,3,4 > repeat-masker.txt
		

		# join requires sorted input
		LANG=en_EN sort -k 1 repeat-masker.txt > repeat-masker-sorted.txt

		# join chromosome identifiers and the data
		LANG=en_EN join -t '	' -1 1 -2 1 chr_map-sorted.txt repeat-masker-sorted.txt > repeat-join.txt

		# remove extra columns
		# FIXME Ensembl uses 1-based coordinates, whereas standard bed file must have 0-based coordinates
		cat repeat-join.txt | cut -d '	' -f 2,3,4 > repeat.bed
	
		# bed to tabix
		cat repeat.bed | sort -k1,1 -k2,2n > repeat-sorted.bed		
		cat repeat-sorted.bed | bgzip > repeat-tabix.bed.gz

		#generate index
                echo $PATH
		tabix -p bed repeat-tabix.bed.gz

		# clean
		rm repeat_feature.txt
		#rm repeat-masker-row.txt repeat-masker-id.txt repeat-sorted.txt repeat-masker-join.txt
		rm repeat-masker.txt  repeat-masker-sorted.txt repeat-join.txt
		rm repeat.bed repeat-sorted.bed


		# Cytoband data


		# sort the actual data
		LANG=en_EN sort -k 2 cytoband-tmp.txt > cytoband-sorted.txt

		# join chromosome identifiers and the data
		LANG=en_EN join -t '	' -1 1 -2 2 chr_map-sorted.txt cytoband-sorted.txt > cytoband-join.txt

		# remove extra columns
		cat cytoband-join.txt | cut -d '	' -f 2,3,4,5,6,7 > cytoband-chr.txt

		# clean
		rm cytoband-sorted.txt cytoband-join.txt "cytoband-tmp.txt"
		rm chr_map-sorted.txt

	#result files: "$2cytoband-chr.txt", "$2repeat-tabix.bed.gz" and "$2repeat-tabix.bed.gz.tbi"
}

chipster_path="0"
#export PATH=${PATH}:/opt/chipster4/comp/modules/admin/shell/:/opt/chipster/tools/emboss/bin/:/opt/chipster/tools/samtools/

ensembl=0
fasta=0
gtf=0
text=0
length=0
version="0.0"
location=$(pwd)
INDEX_BWA=1
INDEX_BOWTIE=1
INDEX_BOWTIE2=1


while [[ $# -ge 1 ]]
do
  case "$1" in
              '-chipster_path')
	      chipster_path="$2"
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
                genome_version="$2"
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
              '-only_bwa')
                INDEX_BOWTIE=0
                INDEX_BOWTIE2=0
                shift
              ;;
	      '-only_bowtie2')
                INDEX_BOWTIE=0
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
tools_path="$chipster_path""/tools"
comp_path="$chipster_path""/comp"
export PATH=${PATH}:$comp_path/modules/admin/shell/:$tools_path/emboss/bin/:$tools_path/samtools/:$tools_path/tabix/tabix-0.2.6/

##
#Retrieve the fasta file
##

#Check if  taxid is used in stead of name

taxid=$(echo $species | tr -d "[a-z,A-Z]" )
species=$(echo $species | sed s/" "/"_"/g )
#test for taxnumber
if [[ "$species" == "$taxid" ]]
then
  species=$(taxget taxon:$taxid -oformat excel -filter | sed s/" "/"_"/g | awk '{print $5}')
  echo "Taxid: $taxid corresponds species: $species"
else
  tax_name=$(echo $species | sed s/"_"/" "/g )
  taxid=$(grep -i "|.$tax_name.|" $tools_path/emboss/share/EMBOSS/data/TAXONOMY/names.dmp  | awk '{print $1}')
fi


echo $species $taxid

#reading the data from ensembl

echo $ensembl

if [[ $ensembl -eq 1 ]]
then
  echo "Retrtieving and indexing genome sequence for $species"

  cd ${tools_path}/genomes/fasta/nochr
  echo ensemblfetch.sh -chipster_path $chipster_path $species
  genome_fasta=$(ensemblfetch.sh -chipster_path $chipster_path $species | tail -1)

  if [[ $genome_fasta == "--------------------------------------------------------------------------------" ]]
  then
    echo Species $species was not found from the Ensembl database.
    exit 1
  fi

  #genome_release=$(echo $genome_fasta | awk -F "." '{print $2}')
  genome_release=$(echo $genome_fasta | awk -F "." '{for (i=2; i<=NF; i++) printf $i"." }' | awk -F ".dna." '{print $1}')

  #Try to find the build name
  wget "http://hgdownload.soe.ucsc.edu/admin/hgcentral.sql"
  # Take only inserts of table dbDB
  cat hgcentral.sql |grep "INSERT INTO \`dbDb\`" > dbDb.txt
  # Take only specific columns
  cat dbDb.txt | cut -d "(" -f 4- |cut -d "'" -f 1,11,16 > release-ucsc.txt
  # Remove some characters
  sed -i 's/)//g' release-ucsc.txt
  sed -i 's/;//g' release-ucsc.txt
  sed -i 's/,//g' release-ucsc.txt
  # Replace slash and quotes with tab
  sed -i 's/\//	/g' release-ucsc.txt #note tab
  sed -i "s/'/	/g" release-ucsc.txt #note tab


  version=$(awk '{if ( $1=="'$genome_release'" ) if ( $NF=="'$taxid'") print $2}' release-ucsc.txt) 
  rm -f release-ucsc.txt
  rm -f hgcentral.sql
  rm -f dbDb.txt
  if [[ $version == "" ]]
  then
    version=$genome_release
  fi

fi



if [[ $fasta -eq 1 ]]
then
  cp $genome_fasta ${tools_path}/genomes/fasta/
  cd ${tools_path}/genomes/fasta
fi

####
#Calculate index file for the fasta file
###
${tools_path}/samtools/samtools faidx $genome_fasta

##
#Check if the fasta file as already been indexed
##

size=$(ls -l $genome_fasta | awk '{print $5}')
checksum=$(md5sum $genome_fasta | awk '{print $1}')

#look for matching size and md5sum

genome_check=$(grep -h $size $tools_path/genomes/genome_list | grep $checksum | awk '{print $1}' | tail -1)


#
if [ ! $genome_check == "" ]; then
  echo "File $genome_fasta has alredy been indexed"
  exit 0
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
  cd ${tools_path}/genomes/gtf
  genome_gtf=$(ensemblfetch.sh -chipster_path $chipster_path -type gtf $species | tail -1 )
  
  process_gtf $genome_gtf

  #FILE_BODY=$(basename $genome_gtf .gtf)
  #cat "$genome_gtf" | cut -f 1,9 --output-delimiter=';' | cut -d ';' -f 1,5  | 	uniq | sed -e 's/; gene_name "/\t/' | sed -e 's/\"//' > "$FILE_BODY.gene.tsv"

fi

if [[ $gtf -eq 1 ]]
then
  if [[ ! -e ${tools_path}/genomes/gtf ]]
  then
     mkdir ${tools_path}/genomes/gtf
  fi  
  cp $location/$genome_gtf ${tools_path}/genomes/gtf/
  FILE_BODY=$(basename $genome_gtf .gtf)
  cat "$location/$genome_gtf" | cut -f 1,9 --output-delimiter=';' | cut -d ';' -f 1,5  | 	uniq | sed -e 's/; gene_name "/	/' | sed -e 's/\"//' > " ${tools_path}/genomes/gtf/$FILE_BODY.gene.tsv"
fi


###
#  get the mysql files
###
if [[ ! -e ${tools_path}/genomes/mysql ]]
then
  mkdir ${tools_path}/genomes/mysql
fi

if [[ $ensembl -eq 1 ]]
then
  cd ${tools_path}/genomes/mysql
  mysql_files=$(ensemblfetch.sh  -chipster_path $chipster_path -type mysql $species | tail -1)
  mysql_dir=$(basename $mysql_files .tar)
  tar xf $mysql_files
  cd $mysql_dir
  gunzip *.gz
  ensembl_mysql .
  #process_mysql_files 
  cd ..
  #ln -s $mysql_dir $species
  rm -f $mysql_files
fi
###


#make bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating BWA indexes for $genome_fasta
  cd $tools_path/bwa_indexes
  if [[ $ensembl -eq 1 ]]
  then
    ln -s ../genomes/fasta/nochr/$genome_fasta $genome_fasta
  else 
    ln -s ../genomes/fasta/$genome_fasta $genome_fasta
  fi
  $tools_path/bwa/bwa index $genome_fasta
else
  echo "Skipping BWA indexing"
fi


#make bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating Bowtie indexes for $genome_fasta 
  cd $tools_path/bowtie/indexes
  if [[ $ensembl -eq 1 ]]
  then
   ln -s ../../genomes/fasta/nochr/$genome_fasta $genome_fasta
  else
   ln -s ../../genomes/fasta/$genome_fasta $genome_fasta
  fi
  bowtie_name=$(basename $genome_fasta .fa)
  $tools_path/bowtie/bowtie-build $genome_fasta $bowtie_name

  #check bowtie2 indexes
  n_genome=$(grep -c "^>" $genome_fasta)
  n_index=$($tools_path/bowtie/bowtie-inspect -n $bowtie_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie indexing of genome $bowtie2_name failed"
     exit 1
  else
    echo Bowtie index of genome $bowtie_name OK
  fi
else
    echo "Skipping Bowtie indexing"
fi





#make bowtie2_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating Bowtie2 indexes for $genome_fasta 

  cd $tools_path/bowtie2/indexes
  if [[ $ensembl -eq 1 ]]
  then
   ln -s ../../genomes/fasta/nochr/$genome_fasta $genome_fasta
  else
    ln -s ../../genomes/fasta/$genome_fasta $genome_fasta
  fi
  bowtie2_name=$(basename $genome_fasta .fa)
  $tools_path/bowtie2/bowtie2-build $genome_fasta $bowtie2_name

  #check bowtie2 indexes
  n_genome=$(grep -c "^>" $genome_fasta)
  n_index=$($tools_path/bowtie2/bowtie2-inspect -n $bowtie2_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie2 indexing of genome $bowtie2_name failed"
     exit 1
  else
    echo Bowtie2 index of genome $bowtie2_name OK
  fi
else
    echo "Skipping Bowtie indexing"
fi


###
#Copy data to genomebowser directory
###

if [[ -e $tools_path/genomebrowser/annotations/$species ]]
then
  echo ""
else
  mkdir $tools_path/genomebrowser/annotations/$species
fi

if [[ -e $tools_path/genomebrowser/annotations/$species/$version ]]
then
  echo ""
else
  mkdir $tools_path/genomebrowser/annotations/$species/$version
fi

echo ""
echo "Copying data to Genome Browser directory:"
echo "$tools_path/genomebrowser/annotations/$species/$version"


cd $tools_path/genomebrowser/annotations/$species/$version

## fasta and fasta index
if [[ $fasta -eq 1 ]]
then
  ln -s ../../../../genomes/fasta/$genome_fasta $genome_fasta
  ln -s ../../../../genomes/fasta/$genome_fasta.fai $genome_fasta.fai
else
  ln -s ../../../../genomes/fasta/nochr/$genome_fasta $genome_fasta
  ln -s ../../../../genomes/fasta/nochr/$genome_fasta.fai $genome_fasta.fai
fi

## gtf
#ln -s  ../../../../genomes/gtf/$genome_gtf 
mv $tools_path/genomes/gtf/$FILE_BODY.gene.tsv ./$FILE_BODY.gene.tsv

##tabix
mv $tools_path/genomes/gtf/$FILE_BODY.tabix.gtf* ./
mv ${tools_path}/genomes/mysql/$mysql_dir/* ./
rm -rf ${tools_path}/genomes/mysql/$mysql_dir
rm -f log repeat-tabix.bed

rm -f karyotype.txt

##genome-v1.yaml

echo species: $species > genome-v1.yaml
echo version: $version >> genome-v1.yaml
echo ensemblBrowserUrl:  >> genome-v1.yaml
echo ucscBrowserUrl: >> genome-v1.yaml
echo sortId: other >> genome-v1.yaml

echo ""
ls -l
echo "---------------------------------------------------------------"



day=$(date)
echo $taxid $species $genome_fasta $version $size $day $checksum >> $tools_path/genomes/genome_list