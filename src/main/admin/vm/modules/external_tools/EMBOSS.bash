##depends:none

source ../installation_files/functions.bash

# EMBOSS, GPL
  apt-get -y install libgd2-noxpm-dev # sudo, emboss needs this to create png images
 	# also vmbin from nic
 	 	 
	##note version in path                                                                                                                                                              
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/EMBOSS-${EMBOSS_VERSION}.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd ${TMPDIR_PATH}/EMBOSS-6.5.7

	# DON'T UNCOMMENT THESE THREE LINES		
  	##wget ftp://emboss.open-bio.org/pub/EMBOSS/fixes/patches/patch-1-11.gz                                                                                                              
 	##gunzip patch-1-11.gz                                                                                                                                                               
	##patch -p1 < patch-1-11                                                                                                                                                             
	
	./configure ${EMBOSS_OPTIONS}
	make
	echo "emboss_update_log.txt" - |make install
	
	# temporarily use prebuilt emboss
	#curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/EMBOSS-${EMBOSS_VERSION}-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
	
	# always create this link
	ln -s EMBOSS-6.5.7 ${TOOLS_PATH}/emboss

# EMBOSS extras

	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/MEME-4.7.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd ${TMPDIR_PATH}/MEME-4.7.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..

# REBASE reference data and indeces	
	cd ${EMBOSS_PATH}/share/EMBOSS/data/REBASE
	#wget_retry ftp://ftp.neb.com/pub/rebase/withrefm.txt
	#wget_retry ftp://ftp.neb.com/pub/rebase/proto.txt
	wget_retry http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/withrefm.txt
	wget_retry http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/proto.txt
	
	../../../../bin/rebaseextract -infile withrefm.txt -protofile proto.txt -equivalences Y
