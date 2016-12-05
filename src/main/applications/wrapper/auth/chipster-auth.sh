#!/bin/bash

# detect chipster home
if [ -h $0 ]
then
    BINARY_REAL_PATH="$(readlink -f "$0")"
else
    BINARY_REAL_PATH="$0"
fi
CHIPSTER_HOME="$(dirname "${BINARY_REAL_PATH}")"

# run
pushd $CHIPSTER_HOME > /dev/null
java -cp '../shared/lib/*' fi.csc.chipster.ChipsterMain authenticator
popd > /dev/null
