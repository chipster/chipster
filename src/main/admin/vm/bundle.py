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
    """
    bundles_yaml = yaml.load(open(filename, "r"))
    logging.debug("load_available_bundles: ")
    pprint(bundles_yaml)
    return bundles_yaml


def load_installed_bundles(filename):
    """
    """
    bundles_yaml = {}
    try:
        bundles_yaml = yaml.load(open(filename, "r"))
        logging.debug("load_installed_bundles: ")
        pprint(bundles_yaml)
    except IOError as e:
        if e.errno == 2:
            print(e)
        else:
            raise
    return bundles_yaml


def save_installed_bundles(filename):
    """
    """
    yaml.dump(installed_bundles, open(filename, "w"), default_flow_style=False)
    logging.debug("save_installed_bundles: Saved!")


def install_bundle(name, version):
    """
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
    """
    if name in installed_bundles:
        return installed_bundles[name]
    else:
        return False


def get_available_bundle(name):
    """
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
    """
    retval = [elem["version"] for elem in get_available_bundle(name)
              if "chipster" in elem and float(elem["chipster"]) <= chipster_version]
    logging.debug("get_compatible_bundle_versions: %s" % retval)
    return retval


def get_available_bundle_version(name, version):
    """
    :type name: str
    :type version: str
    :rtype: list(version)
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

    logging.debug("updated_bundles: %s" % updated_bundles)
    logging.debug("personal_bundles: %s" % personal_bundles)
    logging.debug("deprecated_bundles: %s" % deprecated_bundles)

    return updated_bundles, personal_bundles, deprecated_bundles


def print_available_bundles():
    # TODO: Complete this!
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
            logging.debug(versions)


def create_tree(dst):
    """
    Create tree
    """
    try:
        os.makedirs(os.path.dirname(dst))
    except OSError as e:
        # Tree exists
        if e.errno == 17:
            logging.warning(e)
        else:
            raise
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
        logging.debug("remove_tree(): %s" % dst)
        try:
            os.removedirs(os.path.dirname(dst))
        except OSError as e:
            # Tree doesn't exist
            if e.errno == 2:
                logging.warning(e)
            # Tree not empty
            elif e.errno == 39:
                logging.warning(e)
            else:
                raise
        logging.info("Cleaned tree: %s" % os.path.dirname(dst))

    logging.debug("remove_file(): %s" % dst)
    try:
        os.remove(dst)
        remove_tree(dst)
    except OSError as e:
        # File doesn't exist
        if e.errno == 2:
            logging.warning(e)
        else:
            raise
    logging.info("Removed: %s" % dst)


def transform_bundle(bundle, o_version, n_version):
    """
    Transform an installed bundle version to another version

    Functionality:
        remove = delete files, w/o network traffic needed
        move = move/rename files, w/o network traffic needed
        add = explode package(s) containing new files, w/ network traffic needed
    """

    def get_package_owning_file(tup, bundle, version):
        """
        Get the first package that owns a matching file
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
            shutil.move(refine_path(m[0]), refine_path(m[1]))
        except (OSError, IOError) as e:
            if e.errno == 2:
                logging.warning(e)
            else:
                raise

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
    """
    logging.debug("pkg_name: %s" % pkg_name)
    logging.debug("pkg_values: %s" % pkg_values)

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
        src = file["source"]
        dst = file["destination"]
        checksum = file["checksum"]

        logging.debug("source: %s" % src)
        logging.debug("destination: %s" % dst)
        logging.debug("checksum: %s" % checksum)

        src = tmp_dir + src
        dst = refine_path(dst)

        # Move file into place
        create_tree(dst)
        shutil.move(src, dst)
        logging.info("Moved: %s -> %s" % (src, dst))

    # Loop through symlinks
    if "symlinks" in pkg_values:
        for symlink in pkg_values["symlinks"]:
            create_symlink(symlink["source"], symlink["destination"])

    # Destructively delete temporary directory w/ contents
    shutil.rmtree(tmp_dir)
    logging.info("Temp dir deleted!")
    logging.info("Package %s has exploded!" % pkg_name)


def refine_path(path):
    """
    Refine the path as best as possible

    :type path: str
    :rtype: str
    """
    new_path = path
    if not os.path.isabs(path):
        new_path = installation_path + path
    return new_path


def create_symlink(src, dst):
    logging.debug("source: %s" % src)
    logging.debug("destination: %s" % dst)

    r_dst = refine_path(dst)
    create_tree(r_dst)
    os.symlink(src, r_dst)
    logging.info("Symlinked: %s -> %s" % (refine_path(src), r_dst))


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
                        choices=["install", "uninstall", "update"])  # , "list"])  # ,metavar="action"
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

#    elif args.action == "list":
#        logging.debug("List something!")

#    if name == "all":
#        print_available_bundles()


def diff_bundle(name, version_a, version_b):
    """
    "Calculate" differences between versions of bundle
    NOTE! Should not be dependent on chronology of versions, strictly from version a->b where (a != b)
    What we want to find:
        - added file (checksum in b and not in a)
        - removed file (checksum in a and not in b)
        - moved file (checksum in a and in b, destination in a not equal to that in b)
    """

    def get_checksums(bundle):
        """
        Extract file destination and checksum from bundle contents and return as a set((destination, checksum))
        :type bundle: dict
        :rtype: list(tuple)
        Returns: (source, destination, checksum)
        """
        # logging.debug(bundle)
        checksums = []

        for pkg in bundle["packages"].values():
            for file in pkg["files"]:
                checksums.append((file["source"], file["destination"], file["checksum"]))

        return checksums

    def get_file_for_checksum(checksum, tup):
        """
        """
        return [s[:2] for s in tup if s[2] == checksum][0]

    def get_move_for_checksum(checksum, tuple_a, tuple_b):
        """
        Takes: checksum, (source, destination, checksum), (source, destination, checksum)
        Returns: (old_destination, new_destination), for file with matching checksum
        """
        old = get_file_for_checksum(checksum, tuple_a)[1]
        new = get_file_for_checksum(checksum, tuple_b)[1]

        return old, new

    if float(version_a) == float(version_b):
        raise Exception("This is pointless!")

    content_a = get_available_bundle_version(name, version_a)
    content_b = get_available_bundle_version(name, version_b)
    checksums_a = get_checksums(content_a)
    checksums_b = get_checksums(content_b)
    removed = set([e[2] for e in checksums_a]) - set([e[2] for e in checksums_b])
    added = set([e[2] for e in checksums_b]) - set([e[2] for e in checksums_a])
    diff_a_to_b = set([e[1:] for e in checksums_a]) ^ set([e[1:] for e in checksums_b])
    moved = set([c[1] for c in diff_a_to_b]) - added - removed
    # logging.debug("removed: %s"% removed)
    # logging.debug("added: %s"% added)
    # logging.debug("moved: %s"% moved)

    return ([get_file_for_checksum(a, checksums_b) for a in added],
            [get_file_for_checksum(r, checksums_a) for r in removed],
            [get_move_for_checksum(m, checksums_a, checksums_b) for m in moved])


###########
# Main code
###########

if __name__ == '__main__':
    prog_path = os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
    chipster_version = 2.5
    bundles_file = prog_path + "bundles.yaml"
    installed_file = prog_path + "installed.yaml"
    installation_path = "/tmp/opt/chipster/"
    logging.basicConfig(level=logging.DEBUG)

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
