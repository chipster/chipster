##depends:start/perl_lib.bash

# Install Chipster

cd ${TMPDIR_PATH}/

# use prebuilt version
if [ $CHIPSTER_SOURCE = 'tar' ]; then
  
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/versions/${CHIP_VER}/chipster-${CHIP_VER}.tar.gz | tar -xz

elif [ $CHIPSTER_SOURCE = 'git' ]; then

  #build it from github
  rm -rf chipster
  #git clone https://github.com/chipster/chipster.git
  #cd chipster
  #git checkout $GIT_LABEL 
  curl -s -L http://github.com/chipster/chipster/tarball/$GIT_LABEL/ | tar -xz
  mv chipster-chipster-* chipster
  cd chipster
  
  rm -f keystore.ks
  keytool -genkey -alias chipster -keystore keystore.ks -storepass chipster -keypass chipster -validity 1825 -dname "cn=TEST, ou=TEST, o=TEST, c=TEST"
  echo "chipster" > passfeed; echo "chipster" >> passfeed; 
  ant < passfeed
  cd ..
  mv chipster/dist/chipster-${CHIP_VER}.tar.gz .
  rm -rf chipster
  tar xf chipster-${CHIP_VER}.tar.gz
fi

mv chipster/ ${CHIP_PATH}/

# Make some config "corrections"
sed -i'~' "s:/opt/chipster/:${CHIP_PATH}/:" ${CHIP_PATH}/comp/conf/runtimes.xml
# TODO The below should be made dynamic
sed -i'~' '/<configuration-module moduleId="comp">/a <!-- make compute service access filebroker file repository locally -->\
<entry entryKey="local-filebroker-user-data-path" type="string" \
description="path to local filebrokers user data directory">\
      <value>/opt/chipster/fileserver/file-root/user-data</value>\
</entry>' ${CHIP_PATH}/comp/conf/chipster-config.xml
sed -i'~' "s/#RUN_AS_USER=/RUN_AS_USER=${USERNAME}/" \
    ${CHIP_PATH}/activemq/bin/linux-x86-64/activemq \
    ${CHIP_PATH}/comp/bin/linux-x86-64/chipster-comp \
    ${CHIP_PATH}/auth/bin/linux-x86-64/chipster-auth \
    ${CHIP_PATH}/fileserver/bin/linux-x86-64/chipster-fileserver \
    ${CHIP_PATH}/webstart/bin/linux-x86-64/chipster-webstart \
    ${CHIP_PATH}/manager/bin/linux-x86-64/chipster-manager

# Make update.sh script available
cp ${CHIP_PATH}/admin/vm/update.sh ${CHIP_PATH}/update.sh
chmod u+x ${CHIP_PATH}/update.sh

# If tools foldier doesn't exist, make one
if [ ! -d /mnt/tools ]; then
	mkdir /mnt/tools
fi

# Symlink to tools
ln -s /mnt/tools ${TOOLS_PATH}

# Create user-data and jobs-data symlinks
mkdir ${CHIP_PATH}/fileserver/file-root/
ln -s /scratch/user-data ${CHIP_PATH}/fileserver/file-root/user-data
ln -s /scratch/jobs-data ${CHIP_PATH}/comp/jobs-data

# Symlink to genome browser annotations
mkdir ${CHIP_PATH}/fileserver/file-root/public/
ln -s ${TOOLS_PATH}/genomes/genomebrowser ${CHIP_PATH}/fileserver/file-root/public/annotations

touch ${CHIP_PATH}/auto-config-to-be-run

# TODO Chipster announcements
#mkdir /opt/chipster-announcements
#/opt/chipster-announcements/get-chipster-announcements.sh
# chown -R chipster:chipster /opt/chipster-announcements

#/etc/cron.d/chipster-announcements.crontab
#/etc/update-motd.d/31-vm-instructions
#/etc/update-motd.d/32-announcements
echo "Chipster announcements not done yet"

