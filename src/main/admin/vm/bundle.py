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


def load_available_bundles(filename):
    bundlesYaml = yaml.load(open(filename, "r"))
    print("load_available_bundles: ")
    pprint(bundlesYaml)
    return bundlesYaml

def load_installed_bundles(filename):
    bundlesYaml = {}
    try:
        bundlesYaml = yaml.load(open(filename, "r"))
        print("load_installed_bundles: ")
        pprint(bundlesYaml)
    except IOError as e:
        if e.errno == 2:
            print(e)
        else: raise
    return bundlesYaml

def save_installed_bundles(filename):
    yaml.dump(installed_bundles, open(filename, "w"), default_flow_style=False)
    print("save_installed_bundles: Saved!")

def install_bundle(name, version):
    if is_bundle_installed(name):
        raise Exception("Bundle already installed!")
    # If version not given get the latest compatible
    if not version:
        version = max(get_compatible_bundle_versions(name))
    explode_bundle(name, version)
    installed_bundles[name] = version
    save_installed_bundles(installed_file)

def remove_bundle(name):
    if not is_bundle_installed(name):
        raise Exception("Bundle not installed!")
    # For now version can only be what is already installed
    version = is_bundle_installed(name)
    implode_bundle(name, version)
    del installed_bundles[name]
    save_installed_bundles(installed_file)

def update_bundle(name, n_version):
    # For now o_version can only be what is already installed
    o_version = is_bundle_installed(name)
    if not o_version:
        raise Exception("Bundle not installed!")
    if not n_version:
        n_version = max(get_compatible_bundle_versions(name))
    if n_version == o_version:
        raise Exception("Bundle already this version!")
    remove_bundle(name)
    install_bundle(name, n_version)
    # transform_bundle(name, o_version, n_version)

def is_bundle_installed(name):
    if name in installed_bundles:
        return installed_bundles[name]
    else:
        return False

def get_available_bundle(name):
    retval = None
    if name in available_bundles:
        retval = available_bundles[name]
    print("get_available_bundle:", retval)
    return retval

# def get_available_bundle_versions(name):
#    retval = [elem["version"] for elem in get_available_bundle(name)]
#    print("get_available_bundle_versions:", retval)
#    return retval

def get_compatible_bundle_versions(name):
    retval = [elem["version"] for elem in get_available_bundle(name)
              if "chipster" in elem and float(elem["chipster"]) <= chipster_version]
    print("get_compatible_bundle_versions:", retval)
    return retval

def get_available_bundle_version(name, version):
    retval = [elem for elem in get_available_bundle(name) if elem["version"] == version]
    print("get_available_bundle_version:", retval)
    return retval[0]

def is_bundle_deprecated(name, version):
    retval = [elem["version"] for elem in get_available_bundle(name)
              if "deprecated" in elem and float(elem["deprecated"]) < chipster_version]
    print("is_bundle_deprecated:", retval)
    if retval and retval[0] >= version:
        return True
    else:
        return False

def are_updates_available():
    #
    # Function code
    #
    updated_bundles = {}
    personal_bundles = {}
    deprecated_bundles = {}
    
    if installed_bundles:
        for i_name, i_version in installed_bundles.items():
            print("i_name:", i_name)
            print("i_version:", i_version)
            
            if not get_available_bundle(i_name):
                print("Bundle is personal!\n")
                personal_bundles[i_name] = i_version
            elif is_bundle_deprecated(i_name, i_version):
                print("Bundle is deprecated!\n")
                deprecated_bundles[i_name] = i_version
            elif max(get_compatible_bundle_versions(i_name)) > i_version:
                a_version = max(get_compatible_bundle_versions(i_name))
                print("a_version:", a_version)
                print("Update available!\n")
                updated_bundles[i_name] = a_version
            else:
                print("No update available!\n")
    
    print("updated_bundles:", updated_bundles)
    print("personal_bundles:", personal_bundles)
    print("deprecated_bundles:", deprecated_bundles)
    
    return (updated_bundles, personal_bundles, deprecated_bundles)

# TODO: Complete this!
def print_available_bundles():
    def complement_version_id(tuple):
        print("complement_version_id")
        print("tuple:", tuple)
        if is_bundle_installed(name) == version:
            return version + "*"
        else:
            return version
    
    for bundle, value in available_bundles.items():
        if value:
            versions = sorted([(bundle, version["version"]) for version in value])
            print(versions)

def create_tree(dst):
    '''
    Create tree
    '''
    try: os.makedirs(os.path.dirname(dst))
    except OSError as e:
        # Tree exists
        if e.errno == 17: print(e)
        else: raise
    print("Created tree:", os.path.dirname(dst))

def explode_bundle(name, version):
    '''
    Explode bundle contents into installation path
    '''
    # Blindly assume list can only contain one matching version
    a_version = [b_version for b_version in available_bundles[name] if b_version["version"] == version][0]
    print("version:", a_version)
    
    # Loop through packages
    for pkg_name, pkg_values in a_version["packages"].items():
        explode_package(pkg_name, pkg_values)
    
    print("Bundle", name, version, "has exploded!")
    print()

def implode_bundle(name, version):
    '''
    Implode bundle contents from installation path
    '''
    # Blindly assume list can only contain one matching version
    a_version = [b_version for b_version in available_bundles[name] if b_version["version"] == version][0]
    print("version:", a_version)
    
    # Loop through packages
    for pkg_name, pkg_values in a_version["packages"].items():
        implode_package(pkg_name, pkg_values)
    
    print("Bundle", name, version, "has imploded!")

def remove_tree(dst):
    '''
    Nicely delete only empty directories along path
    '''
    try: os.removedirs(os.path.dirname(dst))
    except OSError as e:
        # Tree doesn't exist
        if e.errno == 2: print(e)
        # Tree not empty
        elif e.errno == 39: print(e)
        else: raise
    print("Cleaned tree:", os.path.dirname(dst))
    print()

def remove_file(dst):
    '''
    Remove file
    '''
    try: os.remove(dst)
    except OSError as e:
        # File doesn't exist
        if e.errno == 2: print(e)
        else: raise
    print("Removed:", dst)

def transform_bundle(bundle, o_version, n_version):
    add, rm, mv = diff_bundle(bundle, o_version, n_version)
    
    print("add", add)
    print("rm", rm)
    print("mv", mv)
    
    # install_file(src, dst)
    for a in add:
        print(a)
    
    for r in rm:
        remove_file(r[1])
        remove_tree(r[1])
    
    for m in mv:
        src = m[0]
        dst = m[1]
        if not os.path.isabs(src):
            src = installation_path + src
        if not os.path.isabs(dst):
            dst = installation_path + dst
        try: shutil.move(src, dst)
        except (OSError, IOError) as e:
            if e.errno == 2: print(e)
            else: raise
    
    print("Bundle", name, o_version, "has transformed into", n_version, "!")

def implode_package(pkg_name, pkg_values):
    print("pkg_name:", pkg_name)
    print("pkg_values:", pkg_values)
    print()
    
    # Loop through files and symlinks
    files = pkg_values["files"]
    if "symlinks" in pkg_values:
        files = files + pkg_values["symlinks"]
    for file in files:
        dst = file["destination"]
        chksum = file["checksum"] if "checksum" in file else None
        
        print("destination:", dst)
        print("checksum:", chksum)
        
        if not os.path.isabs(dst):
            dst = installation_path + dst
        
        remove_file(dst)
        remove_tree(dst)
    print("Package", pkg_name, "has imploded!")

def explode_package(pkg_name, pkg_values):
    print("pkg_name:", pkg_name)
    print("pkg_values:", pkg_values)
    print()
    
    # Download archive
    (f, hm) = urllib.request.urlretrieve(pkg_name)
    print(f)
    print(hm)
    
    # Recognise archive type
    if tarfile.is_tarfile(f):
        print("File is tar (.gz/.bz2)!")
        pf = tarfile.open(f)
    elif zipfile.is_zipfile(f):
        print("File is zip!")
        pf = zipfile.open(f)
    else:
        print("File is unknown!")
        raise Exception("Unknown archive format!")
    
    # Create temporary directory
    tmp_dir = tempfile.mkdtemp() + "/"
    print()
    print("tempdir:", tmp_dir)
    
    # Extract archive
    pf.extractall(tmp_dir)
    pf.close()
    print("Archive unpacked and closed!")
    print()
    
    # Loop through files
    for file in pkg_values["files"]:
        src = file["source"]
        dst = file["destination"]
        chksum = file["checksum"]
        
        print("source:", src)
        print("destination:", dst)
        print("checksum:", chksum)
        
        src = tmp_dir + src
        if not os.path.isabs(dst):
            dst = installation_path + dst
        
        # Move file into place
        create_tree(dst)
        shutil.move(src, dst)
        print("Moved:", src, "->", dst)
        print()
    
    # Loop through symlinks
    if "symlinks" in pkg_values:
        for symlink in pkg_values["symlinks"]:
            src = symlink["source"]
            dst = symlink["destination"]
        
            print("source:", src)
            print("destination:", dst)
            
            if not os.path.isabs(dst):
                dst = installation_path + dst
            
            create_tree(dst)
            os.symlink(src, dst)
            print("Symlinked:", src, "->", dst)
            print()
    
    # Destructively delete temporary directory w/ contents
    shutil.rmtree(tmp_dir)
    print("Temp dir deleted!")
    print("Package", pkg_name, "has exploded!")
    print()

def calculate_checksum(filename):
    '''
    Calculate SHA256 checksum
    '''
    with open(filename, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()

# TODO: Complete this!
def validate_bundle(name, version):
    '''
    Validates the consistency of bundle contents
    '''
    
    # Validate file existens
    os.path.isfile(file)
    
    # Validate symlink existens
    os.path.islink(file)
    
    # Validate checksums
    for file in list_of_files(name):
        if checksum_on_file == calculate_checksum(file):
            print("File is OK!")
        else:
            print("File is corrupted!")

# TODO: Complete this!
def parse_commandline():
    '''
    '''
    
    def get_name_version(str):
        '''
        Return bundle name / version tuple
        '''
        if len(str.split("/")) > 1:
            name, version = str.split("/")
        elif len(str.split("/")) == 1:
            name = str
            version = None
        return (name, version)
    
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
    parser.add_argument("action", type=str, help="Action to perform", choices=["install", "uninstall", "update"])  # , "list"])  # ,metavar="action"
    parser.add_argument("bundle", type=str, help="Bundle <name>/<version> or <keyword>")  # ,metavar="bundle name"
    # parser.add_argument("updates", type=str, help="Check for updates", choices=["check-update"])
    args = parser.parse_args()
    
    name, version = get_name_version(args.bundle)
    print(name, version)
    if args.action == "install":
        print("Install something!")
        install_bundle(name, version)
    elif args.action == "uninstall":
        print("Uninstall something!")
        remove_bundle(name, version)
    elif args.action == "update":
        print("Update something!")
        update_bundle(name, version)
#    elif args.action == "list":
#        print("List something!")
    
#    if name == "all":
#        print_available_bundles()

def diff_bundle(name, version_a, version_b):
    '''
    "Calculate" differences between versions of bundle
    NOTE! Should not be dependent on chronology of versions, strictly from version a->b where (a != b)
    What we want to find:
    - added file (checksum in b and not in a)
    - removed file (checksum in a and not in b)
    - moved file (checksum in a and in b, destination in a not equal to that in b)
    '''
    
    def get_checksums(bundle):
        '''
        Extract file destination and checksum from bundle contents and return as a set((destination, checksum))
        Returns: (source, destination, checksum)
        '''
        # print(bundle)
        checksums = []
        
        for pkg in bundle["packages"].values():
            for file in pkg["files"]:
                checksums.append((file["source"], file["destination"], file["checksum"]))
        
        return checksums
    
    def get_file_for_checksum(checksum, tuple):
        return [s[:2] for s in tuple if s[2] == checksum][0]
    
    def get_move_for_checksum(checksum, tuple_a, tuple_b):
        '''
        Takes: checksum, (source, destination, checksum), (source, destination, checksum)
        Returns: (old_destination, new_destination), for file with matching checksum
        '''
        old = get_file_for_checksum(checksum, tuple_a)[1]
        new = get_file_for_checksum(checksum, tuple_b)[1]
        
        return (old, new)
    
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
    # print("removed:", removed)
    # print("added:", added)
    # print("moved:", moved)
    
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
    
    print("prog_path:", prog_path)
    print("chipster_version:", chipster_version)
    print("bundles_file:", bundles_file)
    print("installed_file:", installed_file)
    print()
    
    available_bundles = load_available_bundles(bundles_file)
    installed_bundles = load_installed_bundles(installed_file)
    
    # update_list, personal_list, deprecate_list = are_updates_available()
    
    # print("calculated checksum:", calculate_checksum("/home/mkarlsso/Downloads/cheatsheet-a4-color.pdf"))
    # explode_bundle("hg19", "1.0")
    # implode_bundle("hg19", "1.0")
    # print(diff_bundle("hg19", "1.0", "1.1"))
    # transform_bundle("hg19", "1.0", "1.1")
    parse_commandline()
