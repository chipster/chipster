##depends:external_tools/EMBOSS.bash

#vienna
                                                                                                                                                                     	curl ftp://emboss.open-bio.org/pub/EMBOSS/VIENNA-1.7.2.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd  ${TMPDIR_PATH}/VIENNA-1.7.2.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..
