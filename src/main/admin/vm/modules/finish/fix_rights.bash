##depends:finish/checksums.bash

## Fix rights 
chown -R -h ${UBUNTU_UID}.${UBUNTU_GID} ${CHIP_PATH}/
if [ -d /mnt/tools/ ]; then
  chown -R ${UBUNTU_UID}.${UBUNTU_GID} /mnt/tools/
  chmod -R go-rwxst,u=rwX,go=rX /mnt/tools/
fi
