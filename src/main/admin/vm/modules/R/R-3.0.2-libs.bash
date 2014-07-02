##depends:R/R-3.0.2.bash

#R libraries, if parallel processing is on, run libs installation with GNU parallel
if [ $parallel == "1" ]; then
  sem --no-notice --tmpdir $paralleldir -j 3 ${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R
else
	${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R
fi

# extra data for zinba R library
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/zinba-extras.tar.gz | tar xz -C ${TOOLS_PATH}

