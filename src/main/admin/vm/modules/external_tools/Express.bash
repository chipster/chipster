##depends:none

#Express, Artistic license 2.0
  #curl http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz | tar -xz -C ${TOOLS_PATH}/
  curl http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/express-1.5.1-linux_x86_64.tgz | tar -xz -C ${TOOLS_PATH}/
  ln -s express-1.5.1-linux_x86_64 ${TOOLS_PATH}/express
