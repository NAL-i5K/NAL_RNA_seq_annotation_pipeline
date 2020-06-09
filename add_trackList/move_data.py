import subprocess
import argparse 
import os
from os import path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-Node1a","--node1_account", help="apollo-nodea account e.g user@apollo-node1,nal.usda.gov",type=str)
parser.add_argument("-bam","--input_bam", help="bam file name",type=str)
parser.add_argument("-s","--Source", help="Source.txt file path",type=str)
args = parser.parse_args()

#get info from Source.txt
with open(args.Source) as f:
    line = f.readline()
    scientific_name = line[:-1]
    line = f.readline()
    assembly_name = line[:-1]
    line = f.readline()
    date = line[:-1]
    Submission = []
    for line in f:
        Submission.append(line[:-1])

gggsss = genus_name[0:3].lower() + species_name[0:3]
temp = scientific_name.split(" ")
genus_name = temp[0]
species_name = temp[1]
bam_file_date = date

#back up trackList on apollo-node1 server
proc_backup = subprocess.Popen(['ssh', args.node1_account, 'cp ' + path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/trackList.json'), path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/trackList.json.bak')])
proc_backup.communicate()

#create new directory on apollo-node1 server
folder_name = genus_name + '-' + species_name + '-RNA-Seq_' + bam_file_date + '_v1.0'
new_dir_path = path.join('/app/data/other_species', gggsss, assembly_name, 'scaffold', 'analyses', folder_name)
proc_mkdir = subprocess.Popen(['ssh', args.node1_account, 'mkdir -p', new_dir_path])
proc_mkdir.communicate()

#rsync files to new directory on apollo-node1 server
proc_rsync = subprocess.Popen(['rsync', '-P', '-r', path.join('/app/data/other_species', gggsss, assembly_name, 'scaffold', 'analyses', folder_name) + '/', args.node1_account + ':' + new_dir_path + '/'])
proc_rsync.communicate()

#scp trackList.json file to apollo-node1 server
proc_scp = subprocess.Popen(['scp', path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/trackList.json'), args.node1_account + ':' + path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/trackList.json')])
proc_scp.communicate()

#create a symlink for analyses folder
proc_ln = subprocess.Popen(['ssh', args.node1_account, 'ln -s', path.join('/app/data/other_species', gggsss, assembly_name, 'scaffold/analyses'), path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data')])
proc_ln.communicate()

#rsync junctions to apollo-node1 server
proc_junction = subprocess.Popen(['rsync', '-P', '-r', path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/tracks', args.input_bam[0:-4] + '_junctions'), args.node1_account + ':' + path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/tracks/')])
proc_junction.communicate()
