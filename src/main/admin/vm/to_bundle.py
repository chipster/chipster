#!/usr/bin/env python3

"""
This tool should create bundle specifications and tar packages from input given on stdin
"""

import logging
import tarfile
import os
import yaml
import bundle


__author__ = 'Mikael Karlsson'


def refine_path(src):
    if not os.path.isabs(src):
        if os.path.exists(chipster_path + src):
            src = chipster_path + src
        elif os.path.exists(tools_path + src):
            src = tools_path + src
        else:
            raise Exception("File path is incomplete!")
    return src


def create_tarball(archive_name, file_list, comp):
    """
    Create a tarball, which is possibly gz or bz2 compressed
    """
    if comp not in ["gz", "bz2"]:
        comp = ""
    tf = tarfile.open(name=archive_name, mode="w:" + comp)
    for file in file_list:
        name_in, name_out = file[0], file[2]
        logging.debug("name_in: %s, name_out: %s" % (name_in, name_out))
        tf.add(name=name_in, arcname=name_out, recursive=False)


def detect_duplicates_and_rename(filename):
    """
    Detect if filename is duplicate and adjust it to allow flat directory structure
    """
    if filename in [fn[2] for fn in file_list]:
        logging.debug("detect_duplicates_and_rename(): filename: %s" % filename)
        if filename[-2] == "." and filename[-1].isdigit():
            filename = "%s%i" % (filename[:-1], int(filename[-1]) + 1)
        else:
            filename += ".1"
    return filename


###########
# Main code
###########

chipster_path = "/tmp/opt/chipster/"
tools_path = chipster_path + "tools/"
chipster_version = 2.5

test_lines = [
    "file1",
    "file2",
    "file3",
    "file4",
    "a/file1"
]
yaml_list = []
yaml_file = "test_spec.yaml"
file_list = []
logging.basicConfig(level=logging.DEBUG)

for file_name in test_lines:
    logging.debug("file name: %s" % file_name)

    # Refine file names to wanted format
    src = refine_path(file_name)
    logging.debug("src: %s" % src)

    dst = refine_path(file_name)
    logging.debug("dst: %s" % dst)

    base_name = os.path.basename(src)
    logging.debug("base name (before): %s" % base_name)

    # Detect duplicates and rename these to "name.[0-9]"
    # This should allow tar archive structure to be flat
    base_name = detect_duplicates_and_rename(base_name)
    logging.debug("base name (after): %s" % base_name)

    # Calculate checksum
    checksum = bundle.calculate_checksum(dst)
    logging.debug("checksum: %s" % checksum)

    # Save all file info to list
    file_list.append((src, dst, base_name, checksum))

logging.debug("list: %s" % file_list)
create_tarball("test.tgz", file_list, "gz")

# Create structure
#abc:
#    - version: x.y
#      chipster: x.y.z
#      deprecated: x.y.z
#      packages:
#        'abc':
#          files:
#            - source: 'abc'
#              destination: 'abc'
#              checksum: '123'
#          symlinks:
#            - source: 'abc/d'
#              destination: 'd'
bundle_name = "abc"
version_number = "1.0"
package_name = bundle_name + ".tgz"
files = [{"source": src, "destination": dst, "checksum": checksum}]
symlinks = [{"source": src, "destination": dst}]
new_packages = {
    package_name: {
        "files": files,
        "symlinks": symlinks
    }
}
new_version = [
    {
        "version": str(version_number),
        "chipster": str(chipster_version),
        "packages": new_packages
    }
]
new_bundle = {bundle_name: new_version}
yaml_list.append(new_bundle)


# Dump yaml
yaml.dump(yaml_list, open(yaml_file, "w"), default_flow_style=False)
