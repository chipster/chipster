##depends:none

# Cufflinks 2, Boost License
cd ${TMPDIR_PATH}/
curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
ln -s cufflinks-2.1.1.Linux_x86_64 ${TOOLS_PATH}/cufflinks2
