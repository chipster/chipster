##depends:none

# primer3, BSD                                                                                                                                                 
	mkdir ${TOOLS_PATH}/primer3
	curl -L http://sourceforge.net/projects/primer3/files/primer3/1.1.4/primer3-1.1.4.tar.gz | tar -xz -C ${TOOLS_PATH}/primer3
	cd ${TOOLS_PATH}/primer3/src/
	make
