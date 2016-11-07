##depends:none                                                                                                                                                                                                 

set -e
#NIC_MIRROR=www.nic.funet.fi
#TMPDIR_PATH=/tmp
#TOOLS_PATH=/opt/chipster/tools
BIOPERL_PATH=${TOOLS_PATH}/bioperl

if [ ! -d $BIOPERL_PATH ]
then
    mkdir /opt/chipster/tools/bioperl
fi


#Berkeley db                                                                                                                                                                                                
mkdir ${BIOPERL_PATH}/berkeleydb
tar xzf /mnt/artefacts/misc/db-6.2.23.tar.gz -C ${BIOPERL_PATH}/berkeleydb/ 
cd ${BIOPERL_PATH}/berkeleydb/db-6.2.23/build_unix/
../dist/configure --prefix=${BIOPERL_PATH}/berkeleydb
make
make install

cd $BIOPERL_PATH

#Perl                                                                                                                                                                                                       
mkdir ${BIOPERL_PATH}/perl
curl -s http://www.cpan.org/src/5.0/perl-5.22.1.tar.gz | tar -xz -C ${TMPDIR_PATH}/
cd ${TMPDIR_PATH}/perl-5.22.1
./Configure -des -Dprefix=${BIOPERL_PATH}/perl/5.22.1
make
make test
make install

OLD_PATH=$PATH
export PATH=${BIOPERL_PATH}/perl/5.22.1/bin:$PATH

#cpanminus                                                                                                                                                                                                  
curl -sL http://search.cpan.org/CPAN/authors/id/M/MI/MIYAGAWA/App-cpanminus-1.7040.tar.gz | tar -xz -C ${BIOPERL_PATH}/perl/
cd ${BIOPERL_PATH}/perl/App-cpanminus-1.7040
perl Makefile.PL
make
make test
make install


#DBI_File can't be installed with cpanm as the Bekrkely DB is not                                                                                                                                           
#installed in default location                                                                                                                                                                              
curl -sL http://search.cpan.org/CPAN/authors/id/P/PM/PMQS/DB_File-1.838.tar.gz | tar -xz -C ${BIOPERL_PATH}/perl/
cd ${BIOPERL_PATH}/perl/DB_File-1.838/
mv config.in config.in.orig
grep -v '^INCLUDE' config.in.orig | grep -v '^LIB' >  config.in
echo "INCLUDE = ${BIOPERL_PATH}/berkeleydb/include" >> config.in
echo "LIB = ${BIOPERL_PATH}/berkeleydb/lib" >> config.in
perl Makefile.PL
make
make test
make install

cd ${BIOPERL_PATH}
cpanm Test::NoWarnings
cpanm CGI
cpanm SVG
cpanm GD --force
cpanm GD::SVG
cpanm XML::Simple
cpanm XML::SAX
cpanm XML::Writer
cpanm Bio::Phylo
#cpanm BioPerl --force
cpanm Bio::Perl --force
cpanm Bio::Graphics --force

#install Virus_Detect                                                                                                                                                                                       
cd /opt/chipster/tools
curl -Ls https://github.com/kentnf/VirusDetect/archive/v1.62.tar.gz | tar -xz -C ${TOOLS_PATH}
cd ${TOOLS_PATH}
ln -s VirusDetect-1.62 virusdetect
cd virusdetect
sed -i -e s/'\/usr\/bin\/perl'/'\/opt\/chipster\/tools\/bioperl\/perl\/5.22.1\/bin\/perl'/g virus_detect.pl
cd bin/
sed -i -e s/'\/usr\/bin\/env perl'/'\/opt\/chipster\/tools\/bioperl\/perl\/5.22.1\/bin\/perl'/g virus_identify.pl

#get databases                                                                                                                                                                                              
cd ${TOOLS_PATH}/virusdetect/databases
db_version=$(curl ftp://bioinfo.bti.cornell.edu/pub/program/VirusDetect/virus_database/ | tail -1 | awk '{ print $NF}' 2> /dev/null)
wget ftp://bioinfo.bti.cornell.edu/pub/program/VirusDetect/virus_database/${db_version}/U95/*
version=$(curl ftp://bioinfo.bti.cornell.edu/pub/program/VirusDetect/virus_database/ | tail -1 | tr -d v | awk '{ print $NF}' 2> /dev/null)
rm plant_*
for type in $(ls *U95.tar.gz | cut -d _ -f1)
do
    tar zxvf ${type}_${version}_U95.tar.gz
    mv ${type}_${version}_U95 $type
    cd $type
    #  for extension in amb ann bwt fai nhr nin nsq pac phr pin psq sa                                                                                                                                      
    for extension in amb ann bwt fai pac sa
    do
        mv *.$extension ../vrl_${type}.$extension
    done
    mv *_prot ../vrl_${type}_prot
    mv vrl_*_U95 ../vrl_${type}
    cd ..
    ../bin/formatdb -i vrl_${type} -p F
    ../bin/formatdb -i vrl_${type}_prot -p T
done

rm -rf ${BIOPERL_PATH}/berkeleydb