#!/usr/bin/env python3

# Let's pretend we're using Python 3
# from __future__ import division, absolute_import, print_function  # , unicode_literals
from pprint import pprint
import argparse
import hashlib
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
import yaml
import zipfile
import logging


def load_available_bundles(filename):
    """
    :type filename: str
    """
    bundles_yaml = yaml.load(open(filename, "r"))
    logging.debug("load_available_bundles: ")
    # pprint(bundles_yaml)
    return bundles_yaml


def load_installed_bundles(filename):
    """
    :type filename: str
    """
    bundles_yaml = {}
    try:
        bundles_yaml = yaml.load(open(filename, "r"))
        logging.debug("load_installed_bundles: ")
        # pprint(bundles_yaml)
    except IOError as e:
        if e.errno == 2:
            print(e)
        else:
            raise
    return bundles_yaml


def save_installed_bundles(filename):
    """
    :type filename: str
    """
    yaml.dump(installed_bundles, open(filename, "w"), default_flow_style=False)
    logging.debug("save_installed_bundles: Saved!")


def install_bundle(name, version):
    """
    :type version: str
    :type name: str
    """
    if is_bundle_installed(name):
        raise Exception("Bundle already installed!")

    # If version not given get the latest compatible
    if not version:
        version = max(get_compatible_bundle_versions(name))
    elif version not in get_compatible_bundle_versions(name):
        raise Exception("Bundle version not compatible or doesn't exist!")
    explode_bundle(name, version)
    installed_bundles[name] = version
    save_installed_bundles(installed_file)


def remove_bundle(name):
    """
    :rtype: None
    :type name: str
    """
    if not is_bundle_installed(name):
        raise Exception("Bundle not installed!")

    # For now version can only be what is already installed
    version = is_bundle_installed(name)
    implode_bundle(name, version)
    del installed_bundles[name]
    save_installed_bundles(installed_file)


def update_bundle(name, n_version):
    """
    :rtype: None
    :type name: str
    :type n_version: str
    """
    # For now o_version can only be what is already installed
    o_version = is_bundle_installed(name)
    if not o_version:
        raise Exception("Bundle not installed!")
    if not n_version:
        n_version = max(get_compatible_bundle_versions(name))
    if n_version == o_version:
        raise Exception("Bundle already this version!")

    # remove_bundle(name)
    # install_bundle(name, n_version)
    transform_bundle(name, o_version, n_version)
    installed_bundles[name] = n_version
    save_installed_bundles(installed_file)


def is_bundle_installed(name):
    """
    :rtype : str or bool
    :type name: str
    """
    if name in installed_bundles:
        return installed_bundles[name]
    else:
        return False


def get_available_bundle(name):
    """
    :rtype : str or None
    :type name: str
    """
    retval = None
    if name in available_bundles:
        retval = available_bundles[name]
    logging.debug("get_available_bundle: %s" % retval)
    return retval

# def get_available_bundle_versions(name):
#   """
#   """
#   retval = [elem["version"] for elem in get_available_bundle(name)]
#   logging.debug("get_available_bundle_versions: %s"% retval)
#   return retval


def get_compatible_bundle_versions(name):
    """
    :rtype: list of str
    :type name: str
    """
    retval = [elem["version"] for elem in get_available_bundle(name)
              if "chipster" in elem and float(elem["chipster"]) <= chipster_version]
    logging.debug("get_compatible_bundle_versions: %s" % retval)
    return retval


def get_available_bundle_version(name, version):
    """
    :type name: str
    :type version: str
    :rtype: dict
    """
    retval = [elem for elem in get_available_bundle(name) if elem["version"] == version]
    logging.debug("get_available_bundle_version: %s" % retval)
    return retval[0]


def is_bundle_deprecated(name, version):
    """
    Check if bundle + version is deprecated

    :type name: str
    :type version: str
    :rtype: bool
    """
    retval = [elem["version"] for elem in get_available_bundle(name)
              if "deprecated" in elem and float(elem["deprecated"]) < chipster_version]
    logging.debug("is_bundle_deprecated: %s" % retval)
    if retval and retval[0] >= version:
        return True
    else:
        return False


def are_updates_available():
    """
    """
    #
    # Function code
    #
    updated_bundles = {}
    personal_bundles = {}
    deprecated_bundles = {}

    if installed_bundles:
        for i_name, i_version in installed_bundles.items():
            logging.debug("i_name: %s" % i_name)
            logging.debug("i_version: %s" % i_version)

            if not get_available_bundle(i_name):
                logging.info("Bundle is personal!")
                personal_bundles[i_name] = i_version
            elif is_bundle_deprecated(i_name, i_version):
                logging.info("Bundle is deprecated!")
                deprecated_bundles[i_name] = i_version
            elif max(get_compatible_bundle_versions(i_name)) > i_version:
                a_version = max(get_compatible_bundle_versions(i_name))
                logging.debug("a_version: %s" % a_version)
                logging.info("Update available!")
                updated_bundles[i_name] = a_version
            else:
                logging.info("No update available!")

    logging.debug("updated_bundles: {}".format(updated_bundles))
    logging.debug("personal_bundles: {}".format(personal_bundles))
    logging.debug("deprecated_bundles: {}".format(deprecated_bundles))

    return updated_bundles, personal_bundles, deprecated_bundles


# TODO: Complete this!
def print_available_bundles():
    """
    """

    def complement_version_id(tup):
        """
        Complement version number given for visual effects
        """
        logging.debug("complement_version_id")
        logging.debug("tuple:", tup)
        if is_bundle_installed(tup[0]) == version:
            return version + "*"
        else:
            return version

    for bundle, value in available_bundles.items():
        if value:
            versions = sorted([(bundle, version["version"]) for version in value])
            print(versions)


def create_tree(dst):
    """
    Create tree
    """
    try:
        os.makedirs(os.path.dirname(dst))
    except OSError as e:
        handle_file_error(e)
    logging.info("Created tree: %s" % os.path.dirname(dst))


def explode_bundle(name, version):
    """
    Explode bundle contents into installation path
    """
    # Blindly assume list can only contain one matching version
    a_version = [b_version for b_version in available_bundles[name] if b_version["version"] == version][0]
    logging.debug("version: %s" % a_version)

    # Loop through packages
    for pkg_name, pkg_values in a_version["packages"].items():
        explode_package(pkg_name, pkg_values)

    logging.info("Bundle %s/%s has exploded!" % (name, version))


def implode_bundle(name, version):
    """
    Implode bundle contents from installation path
    """
    # Blindly assume list can only contain one matching version
    a_version = [b_version for b_version in available_bundles[name] if b_version["version"] == version][0]
    logging.debug("version: %s" % a_version)

    # Loop through packages
    for pkg_name, pkg_values in a_version["packages"].items():
        implode_package(pkg_name, pkg_values)

    logging.info("Bundle %s/%s has imploded!" % (name, version))


def remove_file(dst):
    """
    Remove file
    """

    def remove_tree(dst):
        """
        Nicely delete only empty directories along path
        """
        logging.debug("remove_tree({})".format(dst))
        try:
            os.removedirs(os.path.dirname(dst))
        except OSError as e:
            handle_file_error(e)
        logging.info("Cleaned tree: %s" % os.path.dirname(dst))

    logging.debug("remove_file(): %s" % dst)
    try:
        os.remove(dst)
        remove_tree(dst)
    except OSError as e:
        handle_file_error(e)
    logging.info("Removed: %s" % dst)


def transform_bundle(bundle, o_version, n_version):
    """
    Transform an installed bundle version to another version

    Functionality:
        * remove = delete files, w/o network traffic needed
        * move = move/rename files, w/o network traffic needed
        * add = explode package(s) containing new files, w/ network traffic needed

    :type bundle: str
    :type o_version: str
    :type n_version: str
    """

    def get_package_owning_file(tup, bundle, version):
        """
        Get the first package that owns a matching file
        :type tup: (str,str,str)
        :type bundle: str
        :type version: str
        """
        logging.debug("get_package_owning_file: %s, %s, %s" % (tup, bundle, version))
        for key, values in get_available_bundle_version(bundle, version)["packages"].items():
            for file in values["files"]:
                if file["source"] == tup[0] and file["destination"] == tup[1]:
                    logging.debug("found: %s, %s, %s" % (key, file["source"], file["destination"]))
                    return key, values

    def get_symlinks_for_bundle(name, version):
        """
        Get all symlinks belonging to bundle + version
        :type name: str
        :type version: str
        """
        for x in get_available_bundle_version(name, version)["packages"].values():
            if "symlinks" in x:
                for y in x["symlinks"]:
                    yield y

    add, rm, mv = diff_bundle(bundle, o_version, n_version)

    logging.debug("add %s" % add)
    logging.debug("rm %s" % rm)
    logging.debug("mv %s" % mv)

    for r in rm:
        logging.debug(r)
        remove_file(refine_path(r[1]))

    for m in mv:
        logging.debug(m)
        try:
            shutil.move(refine_path(m[0][1]), refine_path(m[1][1]))
        except (OSError, IOError) as e:
            handle_file_error(e)

    for a in add:
        logging.debug(a)
        pkg_name, pkg_values = get_package_owning_file(a, bundle, n_version)
        logging.debug(pkg_name, pkg_values)
        explode_package(pkg_name, pkg_values)

    # Symlinks, are always removed and added
    [remove_file(refine_path(s["destination"])) for s in get_symlinks_for_bundle(bundle, o_version)]
    [create_symlink(s["source"], refine_path(s["destination"])) for s in get_symlinks_for_bundle(bundle, n_version)]

    logging.info("Bundle %s %s has transformed into %s!" % (bundle, o_version, n_version))


def implode_package(pkg_name, pkg_values):
    """
    :type pkg_name: str
    :type pkg_values: dict
    """
    logging.debug("pkg_name: %s" % pkg_name)
    logging.debug("pkg_values: %s" % pkg_values)

    # Loop through files and symlinks
    files = pkg_values["files"]
    if "symlinks" in pkg_values:
        files = files + pkg_values["symlinks"]
    for file in files:
        dst = file["destination"]
        checksum = file["checksum"] if "checksum" in file else None

        logging.debug("destination: %s" % dst)
        logging.debug("checksum: %s" % checksum)

        remove_file(refine_path(dst))
    logging.info("Package %s has imploded!" % pkg_name)


def explode_package(pkg_name, pkg_values):
    """
    :type pkg_name: str
    :type pkg_values: dict
    """

    def copy_file(src, dst):
        """
        :type src: str
        :type dst: str
        """
        logging.debug("copy_file({})".format(src, dst))

        # Copy file into place
        create_tree(dst)
        # if os.stat(os.path.dirname(src)).st_dev == os.stat(os.path.dirname(dst)).st_dev:
        #     logging.debug("Using link() to copy file!")
        #     os.link(src, dst)
        # else:
        logging.debug("Using copy2() to copy file!")
        shutil.copy2(src, dst)
        # shutil.move(src, dst)
        logging.info("Copied: %s -> %s" % (src, dst))

    logging.debug("explode_package({})".format(pkg_name, pkg_values))

    # Download archive
    (f, hm) = urllib.request.urlretrieve(pkg_name)
    logging.debug(f)
    logging.debug(hm)

    # Recognise archive type
    if tarfile.is_tarfile(f):
        logging.info("File is tar (.gz/.bz2)!")
        pf = tarfile.open(f)
    elif zipfile.is_zipfile(f):
        logging.info("File is zip!")
        pf = zipfile.ZipFile(f)
    else:
        logging.exception("File is unknown!")
        raise Exception("Unknown archive format!")

    # Create temporary directory
    tmp_dir = tempfile.mkdtemp() + "/"
    logging.debug("tempdir: %s" % tmp_dir)

    # Extract archive
    pf.extractall(tmp_dir)
    pf.close()
    logging.info("Archive unpacked and closed!")

    # Loop through files
    for file in pkg_values["files"]:
        copy_file(tmp_dir + file["source"], refine_path(file["destination"]))

    # Loop through symlinks
    if "symlinks" in pkg_values:
        for symlink in pkg_values["symlinks"]:
            create_symlink(symlink["source"], refine_path(symlink["destination"]))

    # Destructively delete temporary directory w/ contents
    shutil.rmtree(tmp_dir)
    logging.info("Temp dir deleted!")
    logging.info("Package %s has exploded!" % pkg_name)


def refine_path(path):
    """
    Refine the path as best as possible
    :type path: str
    """
    new_path = path
    if not os.path.isabs(path):
        new_path = installation_path + path
    return new_path


def create_symlink(src, dst):
    """
    :type src: str
    :type dst: str
    """
    logging.debug("source: %s" % src)
    logging.debug("destination: %s" % dst)

    create_tree(dst)
    try:
        os.symlink(src, dst)
    except OSError as e:
        handle_file_error(e)
    logging.info("Symlinked: %s -> %s" % (dst, src))


def calculate_checksum(filename):
    """
    Calculate SHA256 checksum
    :type filename: str
    :rtype : str
    """
    with open(filename, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def validate_bundle(name, version):
    # TODO: Complete this!
    """
    Validates the consistency of bundle contents
    """

    # Validate file existence
    os.path.isfile(file)

    # Validate symlink existence
    os.path.islink(file)

    # Validate checksums
    for file in list_of_files(name):
        if checksum_on_file == calculate_checksum(file):
            logging.info("File is OK!")
        else:
            logging.warning("File is corrupted!")


def parse_commandline():
    # TODO: Complete this!
    """
    """

    def get_name_version(string):
        """
        Return bundle name / version tuple
        """
        name = None
        version = None
        if len(string.split("/")) > 1:
            name, version = string.split("/")
        elif len(string.split("/")) == 1:
            name = string
        return name, version

    # Actions:
    # install [bundle name/version]
    # upgrade [bundle name/version]
    # uninstall [bundle name]

    # list [installed, available, upgradeable]

    # -v,--verbose

    parser = argparse.ArgumentParser(description="Admin tool for Chipster bundles", epilog="Blah blah blah")
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("-v", "--verbose", action="store_true")
    # group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("action", type=str, help="Action to perform",
                        choices=["install", "uninstall", "update", "list"])  # ,metavar="action"
    parser.add_argument("bundle", type=str, help="Bundle <name>/<version> or <keyword>")  # ,metavar="bundle name"
    # parser.add_argument("updates", type=str, help="Check for updates", choices=["check-update"])
    args = parser.parse_args()

    name, version = get_name_version(args.bundle)
    logging.debug("%s/%s" % (name, version))
    if args.action == "install":
        logging.info("Install something!")
        install_bundle(name, version)
    elif args.action == "uninstall":
        logging.info("Uninstall something!")
        remove_bundle(name)
    elif args.action == "update":
        logging.info("Update something!")
        update_bundle(name, version)
    elif args.action == "list":
        logging.info("List something!")
        print_available_bundles()

#    if name == "all":
#        print_available_bundles()


def diff_bundle(name, version_a, version_b):
    """
    "Calculate" differences between versions of bundle

    NOTE! Should *not* be dependent on chronology of versions, strictly from version 'a'->'b' where ('a' != 'b')

    What we want to find:
        * added file (checksum in 'b' and not in 'a')
        * removed file (checksum in 'a' and not in 'b')
        * moved file (checksum in 'a' and in 'b', destination in 'a' not equal to that in 'b')

    :type name: str
    :type version_a: str
    :type version_b: str
    """

    def get_details(bundle):
        """
        Extract file details from bundle specification and return as a tuple
        :type bundle: dict
        :param bundle: Bundle dictionary
        """
        # logging.debug(bundle)
        for pkg in bundle["packages"].values():
            for file in pkg["files"]:
                yield file["source"], file["destination"], file["checksum"]

    def detect_move(added, removed):
        """
        Detect file move/rename between lists 'added' and 'removed'
        :type added: list of (str,str,str)
        :type removed: list of (str,str,str)
        """

        def old_new(added, removed):
            """
            Sub-function that does the actual work lazily
            :type added: list of (str,str,str)
            :type removed: list of (str,str,str)
            """
            sub_added = list(i[2] for i in added)
            for tup in removed[:]:
                src, dst, checksum = tup
                if checksum in sub_added:
                    i = sub_added.index(checksum)
                    del sub_added[i]
                    removed.remove(tup)
                    yield added.pop(i), tup

        return list(old_new(added, removed)), added, removed

    if float(version_a) == float(version_b):
        raise Exception("This is pointless!")

    spec_a = get_available_bundle_version(name, version_a)
    spec_b = get_available_bundle_version(name, version_b)

    details_a = list(get_details(spec_a))
    details_b = list(get_details(spec_b))
    logging.debug("details_a: %s" % details_a)
    logging.debug("details_b: %s" % details_b)

    moved, added, removed = detect_move(list(b for b in details_b
                                             if b[1:3] not in list(a[1:3] for a in details_a)),
                                        list(a for a in details_a
                                             if a[1:3] not in list(b[1:3] for b in details_b)))
    logging.debug("moved: {}".format(moved))
    logging.debug("added: {}".format(added))
    logging.debug("removed: {}".format(removed))

    return added, removed, moved


def handle_file_error(e):
    """
    :type e: Exception
    :param e: Exception to handle
    """
    # File/Tree doesn't exist
    if e.errno == 2:
        logging.warning(e)
    # File/Tree exists
    elif e.errno == 17:
        logging.warning(e)
    # Tree not empty
    elif e.errno == 39:
        logging.warning(e)
    else:
        raise


###########
# Main code
###########

if __name__ == '__main__':
    prog_path = os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
    chipster_version = 2.6
    bundles_file = prog_path + "bundles.yaml"
    installed_file = prog_path + "installed.yaml"
    installation_path = "/opt/chipster/"
    logging.basicConfig(level=logging.INFO)

    logging.debug("prog_path: %s" % prog_path)
    logging.debug("chipster_version: %s" % chipster_version)
    logging.debug("bundles_file: %s" % bundles_file)
    logging.debug("installed_file: %s" % installed_file)

    available_bundles = load_available_bundles(bundles_file)
    installed_bundles = load_installed_bundles(installed_file)

    # update_list, personal_list, deprecate_list = are_updates_available()

    # logging.debug("calculated checksum: %s"% calculate_checksum("/home/mkarlsso/Downloads/cheatsheet-a4-color.pdf"))
    # explode_bundle("hg19", "1.0")
    # implode_bundle("hg19", "1.0")
    # logging.debug(diff_bundle("hg19", "1.0", "1.1"))
    # transform_bundle("hg19", "1.0", "1.1")
    parse_commandline()
