##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
# sudo aptitude install libhttp-parser-perl
#wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip
wget http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/edirect.zip
unzip -u -q edirect.zip
mv edirect ${TOOLS_PATH}/
rm edirect.zip