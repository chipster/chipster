##depends:none

# see ../docs/VirusDetect-compile.bash

curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/virusdetect/bioperl_ubuntu-16.04.tar.gz | tar -xz -C ${TOOLS_PATH}/
curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/virusdetect/VirusDetect-v1.7_ubuntu-16.04.tar.gz | tar -xz -C ${TOOLS_PATH}/
cd ${TOOLS_PATH}
ln -s VirusDetect-v1.7 virusdetect
