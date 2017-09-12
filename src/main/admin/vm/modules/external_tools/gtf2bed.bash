##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
wget_retry https://raw.githubusercontent.com/timothyjlaurent/GenomicsTools/master/gtf2bed.pl
chmod +x gtf2bed.pl

mkdir -p ${TOOLS_PATH}/gtf2bed
mv gtf2bed.pl ${TOOLS_PATH}/gtf2bed
