##depends:external_tools/EMBOSS.bash

#phylipnew
	curl ftp://emboss.open-bio.org/pub/EMBOSS/PHYLIPNEW-3.69.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd ${TMPDIR_PATH}/PHYLIPNEW-3.69.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
