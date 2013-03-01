#!/usr/bin/env bash

# Fasta files are moved to folder 'real-fasta' and replaced with a symbolink links.
# All annotations files are then compressed to a tar package, which is placed under folder 'annotations'.
#
# The first parameter defines the name of the both packages and the search criteria of the files that belong there. 
# If all necessary files don't start with this package name, the additioanal files can be set with 
# the third parameter.

package () # parameters 1:package 2:fasta 3:files1 4:files2
{
	#Hide the real fasta
	mkdir "real-fasta"
	FASTA=$(basename $2)
	mv "$FASTA" "real-fasta"

	#Replace real file with link to genomes in Chipster virtual machine
	ln -s "../../genomes/fasta/$2" "$FASTA"

	if [ "$3" ] #if threre are additional files
	then
		if [ "$4" ] #if threre are additional files
		then
			tar -cvzf "annotations/$1.tar.gz" $3 $4
		else 
			tar -cvzf "annotations/$1.tar.gz" $3
		fi
	else
		#otherwise, package only files starting with package name
		tar -cvzf "annotations/$1.tar.gz" $1*		
	fi

	#Move original files back so that script setup-genome-browser.sh won't download those again
	mv "real-fasta/$FASTA" .

	rmdir "real-fasta"

}

# report success

exit-trap ()
{
	CODE=$? # keep the relevant exit code 
	if [ ! $CODE -eq "0" ] # there is probably easier way to do this
	then
		echo "Genome Browser annotation packaging FAILED with exit code: $CODE" # how to report problematic command or script line?
	else
		echo "Genome Browser annotation packaging done"
	fi
}


trap exit-trap INT TERM EXIT # execute on exit

echo "Packaging genome browser annotations..."

mkdir annotations

set -e # exit on errors

package "Arabidopsis_lyrata.v.1.0.17" 		"nochr/Arabidopsis_lyrata.v.1.0.16.fa"			"Arabidopsis_lyrata.v.1.0.16.fa"
package "Arabidopsis_thaliana.TAIR10.17" 	"nochr/Arabidopsis_thaliana.TAIR10.17.dna.toplevel.fa"
package "Bos_taurus.UMD3.1.70" 			"nochr/Bos_taurus.UMD3.1.70.dna.toplevel.fa"
package "Canis_familiaris.BROADD2.67" 		"nochr/Canis_familiaris.BROADD2.67.dna.toplevel.fa"
package "Canis_familiaris.CanFam3.1.70" 	"nochr/Canis_familiaris.CanFam3.1.70.dna.toplevel.fa"
package "Gallus_gallus.Gallus_gallus-4.0.pre" 	"nochr/Gallus_gallus.Gallus_gallus-4.0.pre.fa"
package "Gasterosteus_aculeatus.BROADS1.70" 	"nochr/Gasterosteus_aculeatus.BROADS1.70.dna.toplevel.fa"
package "Homo_sapiens.GRCh37.70" 		"nochr/hg19.fa" 					"Homo_sapiens.GRCh37.70*" 	"hg19.fa*"
package "Homo_sapiens.NCBI36.54" 		"nochr/Homo_sapiens.NCBI36.54.dna.toplevel.fa"
package "Human-MT.NC_012920.1"	 		"nochr/Human-MT.NC_012920.1.fa"
package "Mus_musculus.GRCm38.70" 		"nochr/mm10.fa" 					"Mus_musculus.GRCm38.70*" 	"mm10.fa*"
package "Mus_musculus.NCBIM37.67" 		"nochr/mm9.fa" 						"Mus_musculus.NCBIM37.67*" 	"mm9.fa*"
package "N916Ysi" 				"N916Ysi.fa"
package "Ovis_aries_v3.1" 			"ovis_aries_texel.fa" 					"ovis_aries_texel.fa*"
package "Rattus_norvegicus.RGSC3.4.69" 		"nochr/rn4.fa" 						"rn4.fa*"
package "Rattus_norvegicus.Rnor_5.0.70"		"nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa" 
package "R1-RT" 				"R1-RT.fa"
package "Sus_scrofa.Sscrofa10.2.70" 		"nochr/Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa"	"Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa"
package "Vitis_vinifera.IGGP_12x.17" 		"nochr/Vitis_vinifera.IGGP_12x.17.dna.toplevel.fa"

cp contents2.txt annotations/

