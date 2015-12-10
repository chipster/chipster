##depends:R/R-3.0.2.bash


R_VER=3.0.2  
##R libraries, if parallel processing is on, run libs installation with GNU parallel
#if [ $parallel == "1" ]; then
  #  sem --no-notice --tmpdir $paralleldir -j 3 ${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R
#else
	#	${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R
#fi

# extra data for zinba R library, not needed atm
#curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/zinba-extras.tar.gz | tar xz -C ${TOOLS_PATH}

# additions to R-3.0.2-2014-03-03.tar.gz 
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/png-0.1-7-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/QDNAseq-0.99.4-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/
rm -rf ${TOOLS_PATH}/R-3.0.2/lib64/R/library/WECCA
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/WECCA-0.40-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/mgug4852a.db-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/



