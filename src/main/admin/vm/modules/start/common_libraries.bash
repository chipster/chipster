##depends:start/pre_requisites.bash

## Libraries:
# build-essential (only devel)
# gfortran (libgfortran3)
# libcurl4-openssl-dev (libcurl3)
# libglib2.0-dev (libglib2.0-0)
# libglu1-mesa-dev (libglu1-mesa)
# libgsl0-dev (libgsl0ldbl)
# libpng-dev (libpng12-0)
# libreadline-dev (libreadline6) (>libreadline5-dev)
# libxml2-dev (libxml2)
# mesa-common-dev
# tcl-dev (tcl)
# tk-dev (tk)
# xorg-dev (only devel?)
# python-dev (python), for HTSeq
# libnetcdf6, for R
build_tools="yes" # Should tools be built, set to either "yes" or "no"
mode="devel" # Set to either "runtime" or "devel"

if [ $mode == "runtime" ]
then
  ## Runtime:
  aptitude -y --without-recommends install libgfortran3 libcurl3 libglib2.0-0 libglu1-mesa libgsl0ldbl libpng12-0 libreadline6 libxml2 mesa-common-dev tcl tk xorg-dev unixodbc gawk libnetcdf6 cython
elif [ $mode == "devel" ]
then
  ## Devel:
  aptitude -y --without-recommends install build-essential gfortran libcurl4-openssl-dev libglib2.0-dev libglu1-mesa-dev libgsl0-dev libpng-dev libreadline-dev libxml2-dev mesa-common-dev tcl-dev tk-dev xorg-dev python-dev unixodbc-dev libnetcdf-dev openjdk-7-jdk git ant cython
else
  echo "PROBLEM!!"
  exit 1
fi
