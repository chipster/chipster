#!/bin/bash 
# A script to create a new cloud project for the calling user
# This script is called by ../R/create_cloud_project.R
# OH 29.04.2015

#
# EDIT PARAMETERS SPECIFIC FOR LOCAL SERVICE:
#

#source directory to look for files to import into new cloud project
source="/PATH/TO/THE/USERS/DATA"

#chipster command line client credentials. This user needs to have save-as-user rights.
chipster_save_as_user="your-special-user-name-here"
chipster_save_as_user_password="XXXXXXXX"

#
# END OF SPECIFIC PARAMETERS
#

path=`dirname $0`
chipstercli="${path}/chipster-cli.bash"

user=$1
name=$2

filter=""
if [ "$3" != "" ]
then
	filter=`echo "$3" | /bin/sed 's/^filter=//'`
fi
extension=""
if [ "$4" != "" ]
then
	extension=`echo "$4" | /bin/sed 's/^extension=//'`
fi

if [ -f current_time.tmp ]
then
	rm current_time.tmp
fi
date "+%Y%m%d%H%M%S" > current_time.tmp
echo $$ >> current_time.tmp
id=`md5sum current_time.tmp | gawk '{print $1}'`

echo -n "running in "`pwd`$'\r'$'\n'
echo -n "running id is $id"$'\r'$'\n'
echo -n "creating cloud project <$name> for user <$user>"$'\r'$'\n'


if [ ! -d $source ]
then
	echo -n "ERROR: source directory does not exist"$'\r'$'\n'
	exit 1
fi

count=`find $source -maxdepth 1 -name "*${filter}*${extension}" | wc -l | gawk '{print $1}'`
if [ "$count" == "0" ]
then
	echo -n "ERROR: no files matching "'<*'"${filter}"'*'"${extension}"'>'" found in source directory"$'\r'$'\n'
	exit 1
fi

$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose clear-session > /dev/null 2>&1

#for type in fastq bam
#do
	count=`find $source -maxdepth 1 -name "*${filter}*${extension}" | wc -l | gawk '{print $1}'`
	if [ "$count" != "0" ]
	then
		for file in $source/*${filter}*${extension}
		do
			base=`basename $file`
			if [ ! -f $base ]
			then
				echo -n "$file $id"$'\r'$'\n' > $base
				$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose import $base > /dev/null 2>&1
				echo -n "imported $base"$'\r'$'\n'
			else
				echo -n "WARNING: $base already exists, $file not imported"$'\r'$'\n'
			fi
		done
	fi
#done

$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose save-session --cloud $name > /dev/null 2>&1
$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose clear-session > /dev/null 2>&1

#for type in fastq bam
#do
	count=`find $source -maxdepth 1 -name "*${filter}*${extension}" | wc -l | gawk '{print $1}'`
	if [ "$count" != "0" ]
	then
		for file in $source/*${filter}*${extension}
		do
			base=`basename $file`
			count=`find /mnt/data/storage/ -type f -newer current_time.tmp -exec grep -l '^'"$file $id" {} \; | wc -l | gawk '{print $1}'`
			if [ "$count" == "1" ]
			then
				target=`find /mnt/data/storage/ -type f -newer current_time.tmp -exec grep -l '^'"$file $id" {} \;`
				target_base=`basename $target`
				rm $target
				ln -s $file $target
				if [ -f $target.md5 ]
				then
					echo "MAGIC_DO_NOT_CHECK_THIS_FILE_MD5  ${target_base}" > $target.md5
				fi
			else
				echo -n "ERROR: $file in more than one file or not found at all"$'\r'$'\n'
				$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose delete-session --cloud $name > /dev/null 2>&1
				exit 1
			fi
		done
	fi
#done

$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose open-session --cloud $name > /dev/null 2>&1
$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose --save-as-user "$user" save-session --cloud $name > /dev/null 2>&1
$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose delete-session --cloud $name > /dev/null 2>&1
$chipstercli -u "$chipster_save_as_user" -p "$chipster_save_as_user_password" --verbose clear-session > /dev/null 2>&1

echo -n "ready."$'\r'$'\n'$'\r'$'\n'

exit 0 

