##depends:none

source ../installation_files/functions.bash

# primer3, BSD                                                                                                                                                 
	cd ${TMPDIR_PATH}
	mkdir ${TOOLS_PATH}/primer3
	#wget_retry -O primer3-1.1.4.tar.gz http://sourceforge.net/projects/primer3/files/primer3/1.1.4/primer3-1.1.4.tar.gz 
	wget_retry -O primer3-1.1.4.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/primer3-1.1.4.tar.gz 
	tar -xz -C ${TOOLS_PATH}/primer3 -f primer3-1.1.4.tar.gz
	cd ${TOOLS_PATH}/primer3/src/
	make
