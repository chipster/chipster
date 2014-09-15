cd ${TMPDIR_PATH}/
wget http://ftp.gnu.org/gnu/parallel/parallel-20140622.tar.bz2
tar -xvf parallel-20140622.tar.bz2
rm -rf parallel-20140622.tar.bz2
cd parallel-20140622/
./configure
make
make install
cd ..
rm -rf parallel-20140622/

