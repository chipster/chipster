##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
wget https://github.com/enasequence/enaBrowserTools/archive/v1.5.4.tar.gz
tar xf v1.5.4.tar.gz -C ${TOOLS_PATH}/
cd ${TOOLS_PATH}
ln -s enaBrowserTools-1.5.4 enaBrowserTools
