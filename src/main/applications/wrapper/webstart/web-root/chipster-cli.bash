#!/bin/bash

host="http://86.50.168.180:8081"
lib="chipster-current.jar"
conf="chipster-config.xml"
lib_url=$host/lib/$lib
conf_url=$host/$conf

# udpate jar
curl --silent --show-error --remote-name --time-cond $lib $lib_url

# if there are parameters
if (($#))
then
	# set config path
	param="--config $conf_url"
else
	param="-h"
fi


# command for running java
cmd="java -Xmx900M -cp $lib fi.csc.microarray.client.cli.CliClient"

# run
$cmd $param "$@"
