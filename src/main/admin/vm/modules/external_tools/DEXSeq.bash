##depends:none

# DEXSeq
	cd ${TMPDIR_PATH}/
	#	curl -sL http://www.bioconductor.org/packages/release/bioc/src/contrib/DEXSeq_1.8.0.tar.gz | tar -xz
	curl -sL http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.8.0.tar.gz | tar -xz
	mkdir ${TOOLS_PATH}/dexseq-exoncounts
	cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
	cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts
	rm -rf DEXSeq
