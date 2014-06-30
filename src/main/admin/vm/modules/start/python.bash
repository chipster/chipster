##depends:start/pre_requisites.bash

## Python
# !! Anything from PyPI should/shall be installed with pip !!
# python-virtualenv
# virtualenvwrapper
aptitude -y --without-recommends install python-pip
# 2.6
# Needed for MACS?!?, IF REALLY THE CASE IT SHOULD BE RECOMPILED FOR 2.7
aptitude -y install python2.6
