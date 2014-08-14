##depends:none

# EMBOSS, GPL
  apt-get -y install libgd2-noxpm-dev # sudo, emboss needs this to create png images
 	# also vmbin from nic

 	EMBOSS_VERSION=6.5.7
 	EMBOSS_PATH=${TOOLS_PATH}/EMBOSS-${EMBOSS_VERSION}
 	 
	#note version in path                                                                                                                                                              
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/EMBOSS-${EMBOSS_VERSION}.tar.gz | tar -xz -C ${TMPDIR_PATH}/
  cd ${TMPDIR_PATH}/EMBOSS-6.5.7
	
  #wget ftp://emboss.open-bio.org/pub/EMBOSS/fixes/patches/patch-1-11.gz                                                                                                              
 	#gunzip patch-1-11.gz                                                                                                                                                               
	#patch -p1 < patch-1-11                                                                                                                                                             
	
	# used also by other emboss related apps (depend on EMBOSS.bash)
	export EMBOSS_OPTIONS="--prefix=${EMBOSS_PATH}"
	
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	ln -s EMBOSS-6.5.7 ${TOOLS_PATH}/emboss

# EMBOSS extras

	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/MEME-4.7.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd ${TMPDIR_PATH}/MEME-4.7.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..

# REBASE reference data and indeces	
	cd ${EMBOSS_PATH}/share/EMBOSS/data/REBASE
	wget_retry ftp://ftp.neb.com/pub/rebase/withrefm.txt
	wget_retry ftp://ftp.neb.com/pub/rebase/proto.txt
	../../../../bin/rebaseextract -infile withrefm.txt -protofile proto.txt
