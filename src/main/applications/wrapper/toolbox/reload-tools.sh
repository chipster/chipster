#!/bin/bash
echo "triggering tools reload..."

TOOLBOX_HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "toolbox home is ${TOOLBOX_HOME}"

touch ${TOOLBOX_HOME}/.reload/touch-me-to-reload-tools

echo "to check if reload succeeded, tail ${TOOLBOX_HOME}/logs/chipster.log"