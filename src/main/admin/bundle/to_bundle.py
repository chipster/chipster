#!/usr/bin/env python3

"""
This tool should create bundle specifications and tar packages from input given on stdin

NOTES:
    * Create test files e.g with "for n in {1..3}; do dd if=/dev/urandom of=file${n} bs=1M count=1; done"

TODO:
    * Add some form of file duplicate check; warning or symlinking
"""

import argparse
import logging
import pprint
import sys
import tarfile
import os
import yaml
import bundle
import subprocess


__author__ = "Mikael Karlsson <i8myshoes@gmail.com>"


def refine_path(src, tools_path):
    """
    Tries to do some refining and validation on path given
    :rtype: str
    """
    if not os.path.isabs(src):
        if os.path.exists(tools_path + src):
            new_src = tools_path + src
        else:
            raise Exception("File path is incomplete! " + tools_path + src)
    else:
        new_src = src
    return new_src


def create_tarball(archive_name, file_list, compression):
    """
    Create a tarball, which is possibly gz or bz2 compressed
    
    :type archive_name: str
    :type file_list: list
    :type compression: str
    :param archive_name: Archive name
    :param file_list: List of files
    :param compression: Compression schema
        * gz
        * bz2
        * none
    """

    def compress_block(tf):
        logging.debug("compress_block({})".format(tf))
        for file in file_list:
            name_in, name_out = file[0], file[1]
            logging.debug("name_in: %s, name_out: %s" % (name_in, name_out))
            tf.add(name=name_in, arcname=name_out, recursive=False)

    def external_compressor(cmd):
        logging.debug("external_compressor({})".format(cmd))
        try:
            with open(archive_name, "wb") as f, \
                subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=f) as p, \
                tarfile.open(fileobj=p.stdin, mode="w|") as tf:
                compress_block(tf)
        except OSError as e:
            if e.errno == 2:
                logging.warning("Compressor '{}' not found! {}".format(cmd[0], e))
                try:
                    logging.debug("Trying next compressor!")
                    external_compressor(compressors[compression].pop(0))
                except IndexError as e:
                    logging.warning("Out of external compressors! {}".format(e))
                    logging.debug("Switching to internal compressor!")
                    internal_compressor()
            else:
                raise

    def internal_compressor():
        logging.debug("internal_compressor()")
        with tarfile.open(name=archive_name, mode="w:" + compression) as tf:
            compress_block(tf)

    logging.debug("create_tarball({})".format(archive_name, file_list, compression))
    compressors = {
        "gz": [
            ["pigz", "-6"],
            ["gzip", "-6"]
        ],
        "bz2": [
            ["pbzip2", "-6"],
            ["lbzip2", "-6"],
            ["bzip2", "-6"]
        ]
    }

    if compression in compressors:
        external_compressor(compressors[compression].pop(0))
    else:
        internal_compressor()


def detect_duplicates_and_rename(file_name, file_list):
    """
    Detect if filename is duplicate and adjust it to allow flat directory structure
    """
    if file_name in [fn[1] for fn in file_list]:
        logging.debug("detect_duplicates_and_rename(): filename: %s" % file_name)
        if file_name[-2] == "." and file_name[-1].isdigit():
            file_name = "%s%i" % (file_name[:-1], int(file_name[-1]) + 1)
        else:
            file_name += ".1"
    return file_name


def process_file(file_name, file_list):
    """
    Takes a file name and returns a triple (file_name, base_name and checksum)
    """
    logging.debug("file_name: %s" % file_name)
    logging.debug("file_list: %s" % file_list)

    # Get the base name (file name part)
    base_name = os.path.basename(file_name)
    logging.debug("base name (before): %s" % base_name)

    # Detect duplicates and rename these to "name.[0-9]"
    # This should allow tar archive structure to be flat
    base_name = detect_duplicates_and_rename(base_name, file_list)
    logging.debug("base name (after): %s" % base_name)

    # Calculate checksum
    checksum = bundle.calculate_checksum(file_name)
    logging.debug("checksum: %s" % checksum)

    return file_name, base_name, checksum


def create_spec(version, chipster, deprecated, file_list, symlink_list, archive_name):
    """
    This function creates a structure of dicts and lists to be dumped as YAML
    """

    # Files block
    files = []
    for file in file_list:
        files.append(
            {
                "destination": file[0],
                "source": file[1],
                "checksum": file[2]
            }
        )
        logging.debug("files: %s" % files)

    # Symlinks block
    symlinks = []
    for link in symlink_list:
        symlinks.append(
            {
                "source": link[0],
                "destination": link[1]
            }
        )
        logging.debug("symlinks: %s" % symlinks)

    # Package block
    new_packages = {
        archive_name: {
            "files": files,
        }
    }
    if symlinks:
        new_packages[archive_name]["symlinks"] = symlinks

    # Version block
    new_version = {
        "version": str(version),
        "chipster": str(chipster),
        "packages": new_packages
    }
    if deprecated:
        new_version["deprecated"] = str(deprecated)

    return new_version


def main():
    """
    Main function    
    """
    

    yaml_dict = {}
    file_list = []
    symlink_list = []
    logging.basicConfig(level=logging.INFO)

    params = vars(parse_commandline())
    
    if params["tools"]:
        tools_path = params["tools"]
    else:
        tools_path = "/opt/chipster/tools/"
        
    params["archive"] = "{}-{}.tar{}".format(
        params["name"], params["version"],
        "." + params["compression"] if params["compression"] in ("gz", "bz2") else "")
    if not params["file"]:
        params["file"] = "{}-{}.yaml".format(params["name"], params["version"])
    logging.debug("Params: {}".format(params))

    # Process input files
    for file_name in sys.stdin:
        rel_file=file_name.strip()
        abs_file = refine_path(rel_file, tools_path)
        if os.path.islink(abs_file):
            symlink_list.append((os.readlink(abs_file), rel_file))
            logging.debug("symlink_list: %s" % symlink_list)
        elif os.path.isfile(abs_file):
            file_list.append(process_file(file_name=rel_file, file_list=file_list))
            logging.debug("file_list: %s" % file_list)
        elif os.path.isdir(abs_file):
            logging.warning("What are you feeding me!! Directories are rubbish!!")

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

    yaml_dict[params["name"]] = [
        create_spec(version=params["version"],
                    chipster=params["platform"],
                    deprecated=params["deprecated"],
                    file_list=file_list,
                    symlink_list=symlink_list,
                    archive_name=params["prefix"] + params["archive"])
    ]

    pprint.pprint(yaml_dict)
    yaml.dump(yaml_dict, open(params["file"], "w"), default_flow_style=False)
    create_tarball(archive_name=params["archive"],
                   file_list=file_list,
                   compression=params["compression"])


# TODO: Complete this!
def parse_commandline():
    """
    """

    parser = argparse.ArgumentParser(description="Creation tool for Chipster bundles", epilog="Blah blah blah")
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("-v", "--verbose", action="store_true")
    # group.add_argument("-q", "--quiet", action="store_true")
    # parser.add_argument("action",
    #                     type=str,
    #                     help="Action to perform",
    #                     choices=["create"])
    #                     # , "list"])  # ,metavar="action"
    parser.add_argument("-n", "--name",
                        type=str,
                        required=True,
                        help="Bundle <name>")
    # ,metavar="bundle name"
    parser.add_argument("-v", "--version",
                        type=float,
                        required=True,
                        help="Bundle <version>")
    parser.add_argument("-p", "--platform",
                        type=float,
                        required=True,
                        help="Chipster <version>")
    parser.add_argument("-d", "--deprecated",
                        type=float,
                        help="Bundle deprecated since <version>")
    # ,metavar="bundle name"
    parser.add_argument("-c", "--compression",
                        type=str,
                        help="Bundle <compression>",
                        choices=["gz", "bz2", "no"],
                        default="gz")
    parser.add_argument("-f", "--file",
                        type=str,
                        help="Output <file>")
    parser.add_argument("--prefix",
                        type=str,
                        help="package filename prefix, e.g. URL of the bundle repository")
    parser.add_argument("-t", "--tools",
                        type=str,
                        help="Path for converting file names from relative to absolute")
    # ,metavar="bundle name"
    # parser.add_argument("updates", type=str, help="Check for updates", choices=["check-update"])
    args = parser.parse_args()

    return args


###########
# Main code
###########

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt as e:
        logging.warning("Processing interrupted! {}".format(e))
