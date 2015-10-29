##depends:none

#R 3.1.2
R_VER=3.1.2  
#cd ${TMPDIR_PATH}/
#curl -s http://ftp.sunet.se/pub/lang/CRAN/src/base/R-3/R-${R_VER}.tar.gz | tar -xz
#cd R-${R_VER}/
#export MAKEFLAGS=-j
#./configure --prefix=${TOOLS_PATH}/R-${R_VER}
#make
#make install
#echo 'MAKEFLAGS=-j' > ${TOOLS_PATH}/R-${R_VER}/lib/R/etc/Makevars.site # (could also be $HOME/.R/Makevars)
## clean makeflags after R install
#export MAKEFLAGS=
#
#cd ../
#rm -rf R-${R_VER}/

#curl -L http://bio.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.1.2-vmbin/R-3.1.2-2015-02-05.tar.gz | tar -xz -C ${TOOLS_PATH}/
#curl -L http://nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.1.2-vmbin/R-3.1.2-2015-04-08.tar.gz | tar -xz -C ${TOOLS_PATH}/
curl -L http://nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.1.2-vmbin/R-3.1.2-2015-06-01.tar.gz | tar -xz -C ${TOOLS_PATH}/
