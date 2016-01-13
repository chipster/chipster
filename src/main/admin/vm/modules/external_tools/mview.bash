##depends:none

source ../installation_files/functions.bash

# mview

cd ${TMPDIR_PATH}/
#wget http://downloads.sourceforge.net/project/bio-mview/bio-mview/mview-1.59/mview-1.59.tar.bz2
wget http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mview-1.59.tar.bz2
tar xf mview-1.59.tar.bz2 --use-compress-program=pbzip2
cd mview-1.59
sed -e  s/'home\/brown\/HOME\/work\/MView\/dev'/'opt\/chipster\/tools\/mview'/g bin/mview > bin/mview_chipster
mv -f bin/mview_chipster bin/mview
chmod a+x bin/mview
cd ../
mv mview-1.59 ${TOOLS_PATH}/
rm -r mview-1.59.tar.bz2
cd ${TOOLS_PATH}
ln -s mview-1.59 mview
