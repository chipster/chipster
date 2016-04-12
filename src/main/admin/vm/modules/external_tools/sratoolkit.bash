##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
#wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.2-4/sratoolkit.2.4.2-ubuntu64.tar.gz
wget http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/sratoolkit.2.4.2-ubuntu64.tar.gz
cd ${TOOLS_PATH}/
tar zxvf ${TMPDIR_PATH}/sratoolkit.2.4.2-ubuntu64.tar.gz
rm ${TMPDIR_PATH}/sratoolkit.2.4.2-ubuntu64.tar.gz
ln -s sratoolkit.2.4.2-ubuntu64 sratoolkit