##depends:none

#curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.3.2-vmbin/R-3.3.2_ubuntu-16.04_2016-11-28.tar.gz | tar -xz -C ${TOOLS_PATH}/
curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.3.2-vmbin/R-3.3.2_ubuntu-16.04_2017-06-16.tar.lz4 | lz4 -d | tar -x -C ${TOOLS_PATH}/
