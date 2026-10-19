import subprocess
import argparse
import os
import re
from os import path
import shutil

# Rank words skipped when parsing scientific names
RANK_MARKERS = {'subsp.', 'subsp', 'ssp.', 'ssp', 'var.', 'var', 'subvar.',
                'f.', 'forma', 'cf.', 'aff.', 'x', '×'}


def gggsss_type(value):
    """Validate a user-supplied organism code: letters/digits only, stored lowercase."""
    value = value.strip().lower()
    if not re.fullmatch(r'[a-z0-9]+', value):
        raise argparse.ArgumentTypeError(
            'must contain only letters and digits, got "{}"'.format(value))
    return value


def parse_scientific_name(scientific_name):
    """Return [genus, species] or [genus, species, subspecies], cleaned for file names."""
    parts = [p for p in scientific_name.split() if p.lower() not in RANK_MARKERS]
    parts = [re.sub(r'[^A-Za-z0-9]', '', p) for p in parts]
    parts = [p for p in parts if p]
    if len(parts) < 2:
        raise ValueError('Cannot parse genus and species from "{}"'.format(scientific_name))
    if len(parts) > 3:
        parts = parts[:2] + [parts[-1]]   # keep the lowest infraspecific epithet
    return parts


def make_gggsss(name_parts):
    """gggsss for binomials, ggssii (2 letters each) for trinomials; lowercase."""
    if len(name_parts) == 3:
        return ''.join(p[:2] for p in name_parts).lower()
    return (name_parts[0][:3] + name_parts[1][:3]).lower()


parser = argparse.ArgumentParser()
parser.add_argument("-op", "--output_path", help="path of gggsss directory where jbrowse files live (e.g. /app/data/other_species or /project/genomics/username/species_files)", type=str)
parser.add_argument("-j", "--jbrowse_path", help="path of jbrowse container on Scinet", type=str)
parser.add_argument("-bam", "--input_bam", help="bam file name", type=str)
parser.add_argument("-bigwig", "--input_bigwig", help="bigwig file name", type=str)
parser.add_argument("-csi", "--input_csi", help="indexed bam file name, csi format", type=str)
parser.add_argument("-bed", "--input_bed", help="indexed bed file name", type=str)
parser.add_argument("-track", "--input_track", help="trackList.json file path", type=str)
parser.add_argument("-s", "--Source", help="Source.txt file name", type=str)
parser.add_argument("-g", "--gggsss", type=gggsss_type, default=None,
                    help="organism code to use instead of deriving it from the scientific name "
                         "(default: gggsss from genus+species, or ggssii when a subspecies is present). "
                         "Use the same value given to the pipeline's --gggsss.")
parser.add_argument("--overwrite", action="store_true",
                    help="overwrite files that already exist in the analyses folder")
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

#organism code used for both the directory tree and the analyses folder
if args.gggsss:
    gggsss = args.gggsss
    print('Using supplied organism code: {}'.format(gggsss))
else:
    gggsss = make_gggsss(parse_scientific_name(scientific_name))
    print('Derived organism code: {}'.format(gggsss))
folder_name = gggsss + '-RNA-Seq_v1.0'

#warn if the input files were named with a different code
bam = os.path.basename(args.input_bam)
if not bam.startswith(gggsss + '_'):
    print('WARNING: bam file "{}" does not start with "{}_"; check that --gggsss matches the pipeline run'.format(bam, gggsss))

#create the analyses folder (shared by all RNA-seq datasets for this assembly)
new_dir_path = path.join(args.output_path, gggsss, assembly_name, 'scaffold', 'analyses', folder_name)
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
proc_junc = subprocess.Popen(['singularity', 'exec', args.jbrowse_path, 'flatfile-to-json.pl', '--bed', path.join(new_dir_path, bed), '--trackLabel', label + '_junctions', '--config', metadata, '--out', path.join(args.output_path, gggsss, assembly_name, 'jbrowse/data'), '--clientConfig', description])
proc_junc.communicate()
