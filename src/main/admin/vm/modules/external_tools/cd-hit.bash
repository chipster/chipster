##depends:none

source ../installation_files/functions.bash

# cd-hit

cd ${TMPDIR_PATH}/
git clone https://github.com/weizhongli/cdhit
cd cdhit/
make
cd ../
mv cdhit ${TOOLS_PATH}/
