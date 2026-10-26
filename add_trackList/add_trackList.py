import subprocess
import argparse
import os
from os import path
import shutil

from i5k_naming import (species_dir_type, parse_scientific_name, full_name,
                        resolve_species_dir, analyses_folder)

parser = argparse.ArgumentParser()
parser.add_argument("-op", "--output_path", required=True,
                    help="species data root, e.g. /app/data/other_species", type=str)
parser.add_argument("-j", "--jbrowse_path", help="path of jbrowse container on Scinet", type=str)
parser.add_argument("-bam", "--input_bam", help="bam file name", type=str)
parser.add_argument("-bigwig", "--input_bigwig", help="bigwig file name", type=str)
parser.add_argument("-csi", "--input_csi", help="indexed bam file name, csi format", type=str)
parser.add_argument("-bed", "--input_bed", help="indexed bed file name", type=str)
parser.add_argument("-track", "--input_track", help="trackList.json file path", type=str)
parser.add_argument("-s", "--Source", required=True, help="Source.txt file name", type=str)
parser.add_argument("--species_dir", "-g", "--gggsss", dest="species_dir", type=species_dir_type,
                    default=None,
                    help="species directory under --output_path (default: existing legacy gggsss "
                         "directory if present, otherwise Genus_species[_subspecies])")
parser.add_argument("--prefix", type=species_dir_type, default=None,
                    help="file name prefix used by the pipeline (default: Genus_species[_subspecies]). "
                         "Use the same value given to the pipeline's --prefix.")
parser.add_argument("--overwrite", action="store_true",
                    help="overwrite files that already exist in the analyses folder")
parser.add_argument("--dry_run", action="store_true",
                    help="print the resolved paths and exit without copying or editing anything")
args = parser.parse_args()

#get info from Source.txt
with open(args.Source) as f:
    scientific_name = f.readline().rstrip('\n')
    assembly_name = f.readline().rstrip('\n')
    date = f.readline().rstrip('\n')
    Submission = [line.rstrip('\n') for line in f if line.strip()]
print('Scientific name: {}'.format(scientific_name))
print('Assembly name: {}'.format(assembly_name))
print('Annotation files generated date: {}'.format(date))
print('Source: {}'.format(Submission))

#species directory: legacy gggsss if it already exists, otherwise full name
species_dir, reason = resolve_species_dir(args.output_path, scientific_name, args.species_dir)
folder_name = analyses_folder(species_dir)
species_root = path.join(args.output_path, species_dir, assembly_name)
new_dir_path = path.join(species_root, 'scaffold', 'analyses', folder_name)
print('Species directory: {} ({})'.format(species_dir, reason))
print('Analyses folder:   {}'.format(new_dir_path))

#warn if the input files were named with a different prefix
bam = os.path.basename(args.input_bam)
prefix = args.prefix or full_name(parse_scientific_name(scientific_name))
if not bam.startswith(prefix + '_'):
    print('WARNING: bam file "{}" does not start with "{}_"; check that --prefix matches the pipeline run'.format(bam, prefix))

if args.dry_run:
    print('Dry run: nothing copied, trackList.json not modified.')
    raise SystemExit(0)

#create the analyses folder (shared by all RNA-seq datasets for this assembly)
os.makedirs(new_dir_path, exist_ok=True)

#cp files to the analyses folder, refusing to overwrite unless asked
print('copy RNA-annotation output files to scaffold/analyses folder ...')
bigwig = os.path.basename(args.input_bigwig)
csi = os.path.basename(args.input_csi)
bed = os.path.basename(args.input_bed)
copies = [(args.input_bam, bam), (args.input_bigwig, bigwig),
          (args.input_csi, csi), (args.input_bed, bed)]
existing = [name for _, name in copies if path.exists(path.join(new_dir_path, name))]
if existing and not args.overwrite:
    raise SystemExit('These files already exist in {}:\n  {}\nUse --overwrite to replace them.'.format(
        new_dir_path, '\n  '.join(existing)))
for src, name in copies:
    shutil.copyfile(src, path.join(new_dir_path, name))

#backup trackList in apollo directory
shutil.copyfile(args.input_track, args.input_track + '.bak')

#add bam to trackList
print('add bam to trackList ...')
label = os.path.splitext(bam)[0]
data_source_runs = []
for run in Submission:
    data_source_runs.append("<a href='https://www.ncbi.nlm.nih.gov/sra/?term=" + run + "'>" + run + "</a>")
csi_path = path.join('analyses', folder_name, csi)
metadata = '{"csiUrlTemplate":"' + csi_path + '", "unsafePopup": true, "metadata":{"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(data_source_runs) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"Dark red alignments: Mapped portion of read aligned to forward strand<br>Light red alignments: Spliced portion of read aligned to forward strand<br>Dark blue alignments: Mapped portion of read aligned to reverse strand<br>Light blue alignments: Spliced portion of read aligned to reverse strand<br>Red marking - deletion in the read relative to the reference<br>Green marking - insertion in the read relative to the reference<br>Yellow marking  - mismatch (hover over the mismatch to see what the snp is)"}, "type": "WebApollo\/View\/Track\/DraggableAlignments" }'
proc_bam = subprocess.Popen(['singularity', 'exec', args.jbrowse_path, 'add-bam-track.pl', '-i', args.input_track, '-o', args.input_track, '--bam_url', path.join('analyses', folder_name, bam), '--label', label, '--category', 'RNA-Seq/Mapped Reads', '--config', metadata])
proc_bam.communicate()

#add bigwig to trackList
print('add bigwig to trackList ...')
metadata = '{"unsafePopup": true, "metadata":{"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(data_source_runs) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"This track represents an X-Y plot of RNA-Seq coverage."} }'
proc_bw = subprocess.Popen(['singularity', 'exec', args.jbrowse_path, 'add-bw-track.pl', '-i', args.input_track, '-o', args.input_track, '--bw_url', path.join('analyses', folder_name, bigwig), '--label', label + '_coverage', '--plot', '--category', 'RNA-Seq/Coverage Plots', '--pos_color', "#00BFFF", '--config', metadata])
proc_bw.communicate()

#add junction reads to trackList
print('add junction reads  to trackList ...')
metadata = '{"unsafePopup": true, "metadata": {"Analysis provider": "i5k Workspace@NAL",\n"Analysis method": "https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline/",\n"Data source":"' + ','.join(data_source_runs) + '",\n"Publication status":"Analysis: NA; Source data: see individual SRA accessions",\n"Track legend":"Junction reads generated by Hisat2 aligner and regtools. The number of alignments supporting the junction read is listed in blue."}, "category":"RNA-Seq/Splice junctions" }'
description = '{ "description": "score" }'
proc_junc = subprocess.Popen(['singularity', 'exec', args.jbrowse_path, 'flatfile-to-json.pl', '--bed', path.join(new_dir_path, bed), '--trackLabel', label + '_junctions', '--config', metadata, '--out', path.join(species_root, 'jbrowse/data'), '--clientConfig', description])
proc_junc.communicate()
