##depends:none

# Cufflinks 2, Boost License
cd ${TMPDIR_PATH}/
#curl -s http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/

ln -s cufflinks-2.1.1.Linux_x86_64 ${TOOLS_PATH}/cufflinks2
