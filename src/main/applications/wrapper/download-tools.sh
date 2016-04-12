#!/bin/bash                                                                                                     
# Download tools binaries

CHIPSTER_PATH=/opt/chipster
TOOLS_PATH=${CHIPSTER_PATH}/tools
HOST=bio.nic.funet.fi
URL_PATH_PREFIX=/pub/sci/molbio/chipster/dist/virtual_machines

# get Chipster version
CHIPSTER_VERSION=`ls -1 ${CHIPSTER_PATH}/shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`

# construct url 
TOOLS_URL=http://${HOST}/${URL_PATH_PREFIX}/${CHIPSTER_VERSION}/tools/tools.tar.gz

# download
curl ${TOOLS_URL} | tar -xz -C ${TOOLS_PATH}/
