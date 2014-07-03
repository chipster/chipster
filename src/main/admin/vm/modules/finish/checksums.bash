##depends:none

## Create checksums
  cd ${TOOLS_PATH}/
  
  if [ $parallel == "1" ]; then
	parallel --gnu --no-notice -j $jobs -q0 --keep-order -X -N10 $( find . '!' -type d '!' -type l -print0 | xargs -0 sha256sum >> tools.sha256sum )
  else
  	find . '!' -type d '!' -type l -print0 | xargs -0 sha256sum >> tools.sha256sum
fi
