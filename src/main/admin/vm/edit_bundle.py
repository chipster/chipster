#!/usr/bin/env python3

import argparse
import logging
import shutil
import sys
import yaml


def load_yaml(filename):
    """
    :type filename: str
    """
    try:
        bundles_yaml = yaml.load(open(filename, "r"))
    except:
        logging.error("Error in reading yaml file: " + filename)
        raise        
    
    return bundles_yaml

def save_yaml(filename, all_bundles):
    """
    :type filename: str
    """
    yaml.dump(all_bundles, open(filename, "w"), default_flow_style=False)
    logging.debug("save_bundles: Saved!")

def print_bundle(all_bundles, bundle):
    """
    """
    print("%40s" % bundle, end="\t")
    version_list_of_dicts = all_bundles[bundle]                        
    version_list_of_strs = sorted([version["version"] for version in version_list_of_dicts])
    
    for version in version_list_of_strs:                                
        print(version, end="\t")
        
    print("") # new line
                
def list(all_bundles):
    for bundle in all_bundles:            
        print_bundle(all_bundles, bundle)
        
def rename_file_later(old_name, new_name):
    if old_name in files_to_rename:
        raise Exception(
                        "Renaming of file" + old_name + 
                        "to" + new_name + 
                        "failed, because a new name was already stored for it: " + new_name)
    files_to_rename[old_name]=new_name        
    
def apply_file_renames():
    for file in files_to_rename:
        logging.debug("mv from " + file)
        logging.debug("     to " + files_to_rename[file])
        shutil.move(file, files_to_rename[file])
        
def path_to_file(path):
    return path[path.rfind("/") + 1:]

def rename_package(packages, old_name, new_name):
    
    # add a new dict item with the new name and remove the old
    packages[new_name]=packages.pop(old_name)
    rename_file_later(path_to_file(old_name), path_to_file(new_name))
    
def replace_if_unique(replace_in, search_for, replace_with):
    if replace_in.count(search_for)==1:    
        result=replace_in.replace(search_for, replace_with)
        return result
    else:
        raise Exception("Cannot replace " + search_for + 
                        " with " + replace_with + 
                        ", because it isn't unique in " + replace_in)
        
def rename_yaml(yaml_file, search_for, replace_with):
    logging.debug("Search for " + search_for + " in " + yaml_file)    
    if search_for in yaml_file:
        new_yaml = replace_if_unique(yaml_file, search_for, replace_with)
        logging.debug("Rename yaml file " + yaml_file + " to " + new_yaml)                                            
        rename_file_later(yaml_file, new_yaml)      
    
        
def rename_bundle(all_bundles, old_bundle, new_name, yaml_file):
    #rename will modify all_bundles, iterate over a copy of it
    for bundle in all_bundles.copy():  
        for bundle_version in all_bundles[bundle]:
            packages=bundle_version["packages"]
            #rename will modify packages, iterate over a copy of it
            for package in packages.copy():                    
                logging.debug("Rename package " + package)
                new_package = replace_if_unique(package, bundle, new_name)
                rename_package(packages, package, new_package)

        all_bundles[new_name]=all_bundles.pop(bundle)        
        rename_yaml(yaml_file, bundle, new_name)        
                
def set_version(all_bundles, old_bundle, new_version, yaml_file):
      
    for bundle in all_bundles:  
        for bundle_version in all_bundles[bundle]:
                            
            old_version=bundle_version["version"]
            bundle_version["version"]=new_version
                      
            packages = bundle_version["packages"]
            for package in packages.copy():
                new_name = replace_if_unique(package, old_version, new_version)
                rename_package(packages, package, new_name)

            rename_yaml(yaml_file, old_version, new_version)
            
def set_chipster_version(all_bundles, old_bundle, new_version):      
    for bundle in all_bundles:  
        logging.debug("Changing Chipster version for bundle " + bundle)
        for bundle_version in all_bundles[bundle]:
            logging.debug("Chipster version was " + bundle_version["chipster"] + ", will set to " + new_version)            
            bundle_version["chipster"]=new_version
            
def split(all_bundles):
    for bundle in all_bundles:  
        for bundle_version in all_bundles[bundle]:                        
            version=bundle_version["version"]
            split_bundle = {bundle: [bundle_version]}
            filename = bundle + "-" + version + ".yaml"
            save_yaml(filename, split_bundle)
            
def join(files, output):
    
    all_bundles = dict()
    
    for file in files:
        all_bundles.update(load_yaml(file))
        
    save_yaml(output, all_bundles)                        

def process_file(file, args):
    
    print("** Input file: " + file)
        
    all_bundles = load_yaml(file)
    
    selected_bundles = dict();
    unselected_bundles = dict();
    
    if args.bundle:
        for bundle in all_bundles:
            if bundle==args.bundle:  
                selected_bundles[bundle] = all_bundles[bundle]
            else:
                unselected_bundles[bundle] = all_bundles[bundle]
    else:
        selected_bundles = all_bundles
            
         
    if args.split:
        logging.debug("Split something")
        split(selected_bundles)
    elif args.list:
        logging.debug("List something")
        list(selected_bundles)
    elif args.version:
        logging.debug("Set version")
        set_version(selected_bundles, args.bundle, args.version, file)
    elif args.chipster_version:
        logging.debug("Set Chipster version")
        set_chipster_version(selected_bundles, args.bundle, args.chipster_version)
    elif args.name:
        logging.debug("Set name")
        rename_bundle(selected_bundles, args.bundle, args.name, file)                         
        
    if (args.version or args.name or args.chipster):
        save = True        
        if (len(files_to_rename) > 2):
            save = False
            print("\nGoing to rename " + str(len(files_to_rename)) + " files, continue? [y/N]")
            answer = input().lower()
            if answer == "y":
                save = True
            else:
                print("Canceling...")                            
    
        if save:            
            print("Saving changes...")
            all_bundles = selected_bundles
            all_bundles.update(unselected_bundles)            
            save_yaml(file, all_bundles)
            apply_file_renames()
            files_to_rename.clear()
            
def parse_commandline():
    """
    """
        
    end_help="examples:\n" + \
    "  bundles.yaml --list\n" + \
    "  bundles.yaml --split\n" + \
    "  *.yaml --join bundles.yaml\n" + \
    "  bundles.yaml --set_version 0.2\n" + \
    "  bundles.yaml --set_chipster 2.8\n" + \
    "  bundles.yaml --bundle CanFam --set_name Canis_familiaris\n" + \
    "  \n"    

    parser = argparse.ArgumentParser(description="Admin tool for editing Chipster bundles", epilog=end_help, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("input", type=str, nargs="+", help="Yaml file(s) to process")
    
    parser.add_argument("-l", "--list", action="store_true", help="print bundles and versions")
    parser.add_argument("-s", "--split", action="store_true", help="create a separate yaml file for every bundle")
    parser.add_argument("-j", "--join", type=str, dest="output", help="join all input files to single yaml file")    
    parser.add_argument("-v", "--set_version", type=str, dest="version", help="set the provided version number")    
    parser.add_argument("-n", "--set_name", type=str, dest="name", help="rename bundle(s) with the name provided")
    parser.add_argument("-c", "--set_chipster", type=str, dest="chipster_version", help="set required Chipster version")
    parser.add_argument("-b", "--bundle", type=str, help="restrict the operation to single bundle")    
    
    #args = parser.parse_args(["--list", "/home/klemela/tmp/yaml-tool/bundles.yaml"]) # for testing            
    #args = parser.parse_args(["/home/klemela/tmp/yaml-tool/bundles.yaml", "--bundle", "Canis_familiaris.BROADD2","--set_version", "0.2", "-o", "/home/klemela/tmp/yaml-tool/out.yaml"]) # for testing
    #args = parser.parse_args(["/home/klemela/tmp/yaml-tool/Canis_familiaris.BROADD2-0.1.yaml", "--set_name", "Canis_familiaris.NEW_VERSION", "-o", "/home/klemela/tmp/yaml-tool/out.yaml"]) # for testing
    #args = parser.parse_args(["-h"]) # for testing
    args = parser.parse_args() # for command line
        
    if args.output:
        # join requires all files at once
        join(args.input, args.output)
    else:
        # other functions work with the files one by one
        for file in args.input:        
            process_file(file, args)    
        

###########
# Main code
###########

if __name__ == '__main__':
    
    logging.basicConfig(level=logging.INFO)
    
    files_to_rename = dict()    
    parse_commandline()
