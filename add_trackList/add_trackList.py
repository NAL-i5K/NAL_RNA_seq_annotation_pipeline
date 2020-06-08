import subprocess
import argparse
import os
from os import path
import datetime
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-a","--input_account", help="scinet account e.g user@login.scinet.science",type=str)
parser.add_argument("-p","--input_path", help="path of RNA_annotation output files on Scinet",type=str)
parser.add_argument("-bam","--input_bam", help="bam file name",type=str)
parser.add_argument("-bigwig","--input_bigwig", help="bigwig file name",type=str)
parser.add_argument("-bai","--input_bai", help="indexed bam file name",type=str)
parser.add_argument("-bed","--input_bed", help="indexed bed file name",type=str)
parser.add_argument("-track","--input_track", help="trackList.json file path",type=str)
parser.add_argument("-s","--Source", help="Source.txt file name",type=str)
args = parser.parse_args()

#rsync files from Ceres to working directory
proc_rsync = subprocess.Popen(['rsync', '-P', args.input_account + ':' + path.join(args.input_path, args.input_bam) + ' ' + path.join(args.input_path, args.input_bai) + ' ' + path.join(args.input_path, args.input_bigwig) + ' ' + path.join(args.input_path, args.input_bed) + ' ' + path.join(args.input_path, args.Source), '.'])
proc_rsync.communicate()

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
print('Scientific name: {}'.format(scientific_name))
print('Assembly name: {}'.format(assembly_name))
print('Annotation files generated date: {}'.format(date))
print('Source: {}'.format(Submission))

#create new directory
temp = scientific_name.split(" ")
genus_name = temp[0]
species_name = temp[1]
gggsss = genus_name[0:3].lower() + species_name[0:3]
bam_file_date = date
folder_name = genus_name + '-' + species_name + '-RNA-Seq_' + bam_file_date + '_v1.0'
new_dir_path = path.join('/app/data/other_species', gggsss, assembly_name, 'scaffold', 'analyses', folder_name)
os.makedirs(new_dir_path)

#cp files to new directory
print('copy RNA-annotation output files to scaffold/analyses folder ...')
shutil.copyfile(args.input_bam, path.join(new_dir_path, args.input_bam))
shutil.copyfile(args.input_bigwig, path.join(new_dir_path, args.input_bigwig))
shutil.copyfile(args.input_bai, path.join(new_dir_path, args.input_bai))
shutil.copyfile(args.input_bed, path.join(new_dir_path, args.input_bed))

#create symlink to bam, indexed and bigwig files
src_dir = path.join('/app/data/other_species', gggsss, assembly_name, 'scaffold', 'analyses')
dst_dir = path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data/analyses')
if not path.exists(dst_dir):
    print('Create a symlink')
    os.symlink(src_dir, dst_dir)

#backup trackList in apollo directory
shutil.copyfile(path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data', 'trackList.json'), path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data', 'trackList.json.bak'))

#add bam to trackList
print('add bam to trackList ...')
label = os.path.splitext(os.path.basename(args.input_bam))[0]
metadata = '{"metadata":{"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(Submission) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"Dark red alignments: Mapped portion of read aligned to forward strand<br>Light red alignments: Spliced portion of read aligned to forward strand<br>Dark blue alignments: Mapped portion of read aligned to reverse strand<br>Light blue alignments: Spliced portion of read aligned to reverse strand<br>Red marking - deletion in the read relative to the reference<br>Green marking - insertion in the read relative to the reference<br>Yellow marking  - mismatch (hover over the mismatch to see what the snp is)"}, "type": "WebApollo\/View\/Track\/DraggableAlignments" }'
proc_bam = subprocess.Popen(['add-bam-track.pl', '-i', args.input_track, '-o', args.input_track, '--bam_url', path.join('analyses', folder_name, args.input_bam), '--label', label[0:-10] + bam_file_date, '--category', 'RNA-Seq/Mapped Reads', '--config', metadata])
proc_bam.communicate()

#add bigwig to trackList
print('add bigwig to trackList ...')
metadata = '{"metadata":{"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(Submission) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"This track represents an X-Y plot of RNA-Seq coverage."} }'
proc_bw = subprocess.Popen(['add-bw-track.pl', '-i', args.input_track, '-o', args.input_track, '--bw_url', path.join('analyses', folder_name, args.input_bigwig), '--label', label[0:-10] + bam_file_date + '_coverage', '--plot', '--category', 'RNA-Seq/Coverage Plots', '--pos_color', "#00BFFF", '--config', metadata])
proc_bw.communicate()

#add juntion reads to trackList
print('add junction reads  to trackList ...')
metadata = '{"style":{"showLabels": false}, "metadata": {"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(Submission) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"Junction reads generated by Hisat2 aligner and regtools"}, "category":"RNA-Seq/Splice junctions" }'
proc_junc = subprocess.Popen(['flatfile-to-json.pl', '--bed', path.join(new_dir_path, args.input_bed), '--trackLabel', label[0:-10] + bam_file_date + '_junctions', '--config', metadata, '--out', path.join('/app/data/other_species', gggsss, assembly_name, 'jbrowse/data')])
proc_junc.communicate()
