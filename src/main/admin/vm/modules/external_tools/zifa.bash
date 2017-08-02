##depends:git, python

# ZIFA
cd ${TMPDIR_PATH}/
git clone https://github.com/epierson9/ZIFA
cd ZIFA
# TODO: on commandline this requires sudo
${TOOLS_PATH}/Python-2.7.12/bin/python setup.py install
cd ..
rm -rf ZIFA
