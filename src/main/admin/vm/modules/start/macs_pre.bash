##depends:start/chipster_main.bash

# MACS, Artistic license
# part 1
cd ${TMPDIR_PATH}/
curl -s https://cloud.github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz | tar -xz
cd MACS-1.4.2
python setup.py install
cd ..
rm -rf MACS-1.4.2
rm -f MACS-1.4.2-1.tar.gz
