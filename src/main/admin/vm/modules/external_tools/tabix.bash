##depends:none

# tabix
  # curl -sL http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download | tar xj -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/tabix-0.2.6-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s tabix-0.2.6 ${TOOLS_PATH}/tabix 
