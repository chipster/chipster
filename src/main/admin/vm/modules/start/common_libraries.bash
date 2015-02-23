##depends:start/pre_requisites.bash

## Libraries:
# python-dev (python), for HTSeq
# libnetcdf6, for R

aptitude -y --without-recommends install build-essential gfortran libcurl4-openssl-dev libglib2.0-dev libglu1-mesa-dev libgsl0-dev libpng-dev libreadline-dev libxml2-dev mesa-common-dev tcl-dev tk-dev xorg-dev python-dev unixodbc-dev libnetcdf-dev openjdk-7-jdk git ant cython
