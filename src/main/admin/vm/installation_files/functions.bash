#!/usr/bin/env bash

# Installation functions

# Function to install one module
# It uses flagfiles to check if module is already installed

install_module()
{

	flagfile="../installation_files/flags/$1"
	
	# Make checksum for module file
	checksum=$( sha256sum $1 )
	
	# If flagfile exists, take checksum from the file
if [ -f $flagfile ]; then
	flagsum=$( head -n 1 $flagfile )
else 
	
	# If not, set flagsum to "None" to prevent unbound variable

	flagsum="NONE"		
fi

# If force is used, don't check flagfiles
if [[ "$force" == "0" ]]; then

		
		# Check if flagfiles checksum matches
		# don't install module
		if [ "$flagsum" == "$checksum" ]; then
			echo "$1 already installed"
		else 

		# Install module
		(
		
			bash -veu $1
		
		)
	
		#If installed correctly, make sha256sum & flagfile
		sha256sum $1 > $flagfile
	fi
else
	(
    	# Install module
		echo "** Install $1"
		if bash -eu $1; then
			echo "** Install $1 done"
		else
			echo "** Install $1 FAILED"
			exit 1
		fi
	)
	
	#If installed correctly, make sha256sum & flagfile
	sha256sum $1 > $flagfile
	
fi
	
}

export -f install_module

# Function to call install_module for all lines in the file

install_file() {
	if [ $parallel == "1" ]; then
		parallel --gnu --no-notice -v -j $jobs -a $1_parallel install_module
	fi		
	while read line; do
		install_module "$line"
	done < $1
		
}

export -f install_file


# Function to do a hard retry with wget
	
wget_retry()
{
	max_retries=10
	i=1
	while ! wget -t 1 -T 30 $*; do
		if [ $i -ge $max_retries ]; then
			echo "maximum number of retries exceeded"
			return 1
		fi
		sleep 5
		i=$((i+1))
		echo "wget failed, retrying ($i/$max_retries)"
		
	done
}

export -f wget_retry

# Function to install modules which depedency is installed module

install_dep() {
	
	directory=$( echo $1 | cut -d "/" -f1 )
	modules=( $(find $directory/*.bash) )
	
	for i in "${modules[@]}"; do
		
		# take one module depency
		depency=( $(cut -d : -f2 $i) )
		
		# If module depends from $1, install it
		if [[ "$depency" == "$1" ]]; then
			cd ..
			cd ..
			bash install-module.sh $i
			
		fi
	done
}

export -f install_dep

# Function to create flag-directories
create_flag_dirs() {
	while read line; do
		mkdir -p installation_files/flags/$line
	done < $1
}

export -f create_flag_dirs

read_dirs() {
	while read line; do
		folders+=( $line )
	done < $1
	export folders
}

export -f read_dirs

set -e

	