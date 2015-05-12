# TOOL extract_gzip_tar.py: "Extract .tar.gz file" (Extract a gzipped tar file, which usually has a file extension .tar.gz)
# INPUT input_file: ".tar.gz file" TYPE GENERIC (Gzip compressed tar file)
# OUTPUT output_file{...}: "Extracted file(s)"

import tarfile
import posixpath
from tool_utils import *

def main():

    input_names = read_input_definitions()
    output_names = {}

    # benchmark this in comparison to gunzip and pigz
    infile = tarfile.open('input_file', 'r')

    for member in infile.getmembers():
        if not member.isfile():
            print 'skipping, not a file: ' + member.name
            continue
        # fixed names for output files
        output_file = 'output_file' + str(len(output_names) + 1)
        # remove paths from dataset names, because those aren't supported in client
        output_names[output_file] = posixpath.basename(member.name)
        # extract without paths, because those could point to parent directories
        member.name = output_file
        infile.extract(member, '.')

    infile.close()

    # set dataset names
    write_output_definitions(output_names)

if __name__ == "__main__":
    main()
