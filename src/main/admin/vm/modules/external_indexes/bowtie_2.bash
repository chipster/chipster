##depends:none

# Bowtie2 indexes, built for Chipster
  cd ${TMPDIR_PATH}/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/	
  	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
	
ln -s -t ${TOOLS_PATH}/bowtie2/indexes ../../genomes/fasta/nochr
rm -f ${TOOLS_PATH}/bowtie2/indexes/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa
ln -s ../../genomes/fasta/nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa ${TOOLS_PATH}/bowtie2/indexes

# Tophat indexes
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/hg19.ti.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/

