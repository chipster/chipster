##depends:start/macs2_pre.bash

# HTSeq, GPL v3 or later
# part 1
source ../installation_files/functions.bash
pip install HTSeq==0.6.1
wget_retry -O /usr/local/bin/htseq-count_chr http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
chmod 755 /usr/local/bin/htseq-count_chr
wget_retry -O /usr/local/lib/python2.7/dist-packages/HTSeq/scripts/count_chr.py http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr_v2.py
