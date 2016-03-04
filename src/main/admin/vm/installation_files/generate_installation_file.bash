# Bash script to list&sort files to installation file
# takes two arguments: foldier and wanted name of the installation file


# Arrays for sorted depencies and temporary depency_array
sort_array=()
depency_array=()
parallel_array=()


#Array function:

arraySort() {
	
	# save given modules to array
	array=( "$@" )
	
	# for the safety, unset depency array
	unset depency_array

		 
# for array size
	for i in "${array[@]}"; do
	
	# debug
	echo "$i"
		
	# take one module depency
	depency=( $(cut -d : -f2 $i) )
	
	
	# remove taken module from array
	array=( "${array[@]/$i}" )
	
	
	# if depency is none, put first to sort_array
	if [[ $depency == "none" ]]; then
		sort_array+=( $i )
		if [ $parallel == "1" ]; then
			parallel_array+=( $i )
		fi
	
	# if there is depency
	else
		
		
		#Help variable
		dep_not_found=1
		
		# look if depency is found from sort_array
		for dep in "${sort_array[@]}"; do
			
    		if [[ "$dep" == "$depency" ]]; then
    		
			# put module to sort_array's end
			sort_array+=( $i )
			
			dep_not_found=0
    		fi
		done
		
		# if not, check if depency exist and put module to depency_array
		if [[ "$dep_not_found" -eq 1 ]]; then
			# Add module to depency_array
			depency_array+=( $i )
		fi
		
	fi
	done
		
	#if depency_array != empty
	if [ ${#depency_array[@]} -ne 0 ]; then
	
		#recurse Array function with depency_array
		arraySort "${depency_array[@]}"
	fi
}

# Move to working directory
cd modules/

# List files from directory to the array
modules=( $(find $1*.bash) )


# Call arraySort function with modules
arraySort "${modules[@]}"


# Save filenames from array to $1-given file
# if parallel in use, do two files
	if [ $parallel == "1" ]; then
		for i in "${parallel_array[@]}"; do
			sort_array=( ${sort_array[@]/$i} )
		done
		printf "%s\n" "${parallel_array[@]}" > $2_parallel
	fi

printf "%s\n" "${sort_array[@]}" > $2
sed '/^$/d' $2 > $2.out
mv  $2.out $2

# Return to home foldier
cd ..
