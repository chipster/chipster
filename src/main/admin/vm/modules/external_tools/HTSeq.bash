##depends:none

# HTSeq, GPL v3 or later
  # part 2
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/htseq/
  ln -s /usr/local/bin/htseq-qa ${TOOLS_PATH}/htseq/htseq-qa
  ln -s /usr/local/bin/htseq-count ${TOOLS_PATH}/htseq/htseq-count
  ln -s /usr/local/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr
