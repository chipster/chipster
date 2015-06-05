# TOOL gunzip.py: "Extract .gz file" (Extract a gzip file, which usually has a file extension .gz)
# INPUT input_file: "Gzip file" TYPE GENERIC (Gzip compressed file)
# OUTPUT output_file: "Extracted file"

import gzip
import shutil
from tool_utils import *

def main():

    infile = gzip.open('input_file', 'rb')

    with open('output_file', 'wb') as outfile:
        # copy in chunks
        shutil.copyfileobj(infile, outfile)
    infile.close()

    # set dataset name
    input_name = read_input_definitions()['input_file']
    output_names = {'output_file': remove_postfix(input_name, '.gz')}

    write_output_definitions(output_names)

if __name__ == "__main__":
    main()
