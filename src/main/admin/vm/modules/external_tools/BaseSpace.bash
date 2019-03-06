##depends:none

source ../installation_files/functions.bash

mkdir  ${TOOLS_PATH}/BaseSpace_cli-0.10.7
mkdir  ${TOOLS_PATH}/BaseSpace_cli-0.10.7/bin

wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O  ${TOOLS_PATH}/BaseSpace_cli-0.10.7/bin/bs
wget https://api.bintray.com/content/basespace/BaseSpace-Copy-BIN/\$latest/linux/bscp?bt_package=develop -O  ${TOOLS_PATH}/BaseSpace_cli-0.10.7/bin/bs-cp

cd ${TOOLS_PATH}
ln -s BaseSpace_cli-0.10.7 bs

