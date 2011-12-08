#!/bin/bash                                                                                                     
if [ "$#" -eq 0 ]
then
    java -jar shared/lib/chipster-config-tool.jar configure
elif [ "$#" -eq 1 ]
then
    if [ "$1" = "auto" ]
    then
        java -jar shared/lib/chipster-config-tool.jar auto-configure
    else
        java -jar shared/lib/chipster-config-tool.jar simple-configure $1
    fi
else
    echo "Usage: ./configure.sh <hostname or address> | auto"
fi


