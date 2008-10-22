#!/bin/sh

INSTALL_SCRIPT="shared/lib/install.R"
R_PARAMS="--vanilla"
DEFAULT_R_COMMAND="R"
R_VERSION="2.6.1"

echo "Chipster R package installer"
echo
echo "No changes are written before you verify them"
echo

# CHECK

echo "Before proceeding you must have:"
echo "  1) R version $R_VERSION installed" 
echo "  2) permissions for installing R packages (usually you have to be root)"
echo "Are these conditions met [yes/no]?"
read answer

if [ "$answer" != "yes" ]; then
    echo "User decided to abort, no changes made"
    exit
fi

# GATHER INFO

echo "Please specify R command [R]: "
read command

command=${command## } # trim (remove leading white space)
command=${command%% } # trim (remove trailing white space)
if [ "$command" = "" ]; then
    command=$DEFAULT_R_COMMAND
fi
echo


# VERIFY CHANGES

echo "Using R command \"$command\""
echo "Using the following install script ($INSTALL_SCRIPT):"
cat $INSTALL_SCRIPT
echo 

echo "Please verify the install script. Should the script be run [yes/no]?"
read answer

if [ "$answer" != "yes" ]; then
    echo "User decided to abort, no changes made"
    exit
fi
       

# MAKE CHANGES

cat $INSTALL_SCRIPT | $command $R_PARAMS
echo "All changes successfully made!"
