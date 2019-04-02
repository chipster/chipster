##depends:none

source ../installation_files/functions.bash

#https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
#curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/sratoolkit.current-ubuntu64.tar.gz | tar -xz -C ${TOOLS_PATH}/
curl -L http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.4/sratoolkit.2.9.4-ubuntu64.tar.gz | tar -xz -C ${TOOLS_PATH}/
cd ${TOOLS_PATH}
ln -s sratoolkit.2.9.4-ubuntu64 sratoolkit
