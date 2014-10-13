##depends:start/macs_pre.bash

# MACS2, BSD licence
# part 1

cd $TMPDIR_PATH
#curl -s 'https://cloud.github.com/downloads/taoliu/MACS/MACS-2.0.9-1.tar.gz' | tar -xz
curl -s 'https://pypi.python.org/packages/source/M/MACS2/MACS2-2.1.0.20140616.tar.gz' | tar -xz

cd MACS2-2.1.0.20140616/
sudo python setup.py install
cd $TMPDIR_PATH
rm -rf MACS2-2.1.0.20140616
