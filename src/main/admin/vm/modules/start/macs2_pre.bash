##depends:start/macs_pre.bash

# MACS2, BSD licence
# part 1

cd $TMPDIR_PATH
curl -s 'https://cloud.github.com/downloads/taoliu/MACS/MACS-2.0.9-1.tar.gz' | tar -xz
cd MACS-2.0.9/
sudo python setup.py install
cd ..
rm -rf MACS-2.0.9
