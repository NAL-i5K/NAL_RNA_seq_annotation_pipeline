import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxid", help="your tax_id", type=str)
parser.add_argument(
    "-o",
    "--output",
    help=("output file name, ",
          "if not specified, will generate output tsv file ",
          "in current folder with name [tax_id].tsv"),
    type=str)
args = parser.parse_args()
tax_id = args.taxid
print('Processing tax id:' + args.taxid)

proc_esearch = subprocess.Popen(
    ['esearch', '-db', 'taxonomy', '-query', tax_id + '[uid]'],
    stdout=subprocess.PIPE)
proc_elink = subprocess.Popen(
    ['elink', '-target', 'sra'],
    stdin=proc_esearch.stdout,
    stdout=subprocess.PIPE)
proc_efilter = subprocess.Popen(
    ['efilter', '-query', 'rna seq[stra] AND Transcriptomic[src]'],
    stdin=proc_elink.stdout,
    stdout=subprocess.PIPE)
proc_efetch = subprocess.Popen(
    ['efetch', '-format', 'runinfo', '-mode', 'xml'],
    stdin=proc_efilter.stdout,
    stdout=subprocess.PIPE)
proc_xtract = subprocess.Popen(
    [
        'xtract', '-pattern', 'Row', '-def', 'N/A', '-element', 'Run',
        'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates',
        'avgLength', 'size_MB', 'download_path', 'Experiment', 'LibraryName',
        'LibraryStrategy', 'LibrarySelection', 'LibrarySource',
        'LibraryLayout', 'InsertSize', 'InsertDev', 'Platform', 'Model',
        'SRAStudy', 'BioProject', 'Study_Pubmed_id', 'ProjectID', 'Sample',
        'BioSample', 'SampleType', 'TaxID', 'ScientificName', 'SampleName',
        'Sex', 'Tumor', 'Submission', 'Consent', 'RunHash', 'ReadHash'
    ],
    stdin=proc_efetch.stdout,
    stdout=subprocess.PIPE)

output_header = '\t'.join([
    'Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates',
    'avgLength', 'size_MB', 'download_path', 'Experiment', 'LibraryName',
    'LibraryStrategy', 'LibrarySelection', 'LibrarySource', 'LibraryLayout',
    'InsertSize', 'InsertDev', 'Platform', 'Model', 'SRAStudy', 'BioProject',
    'Study_Pubmed_id', 'ProjectID', 'Sample', 'BioSample', 'SampleType',
    'TaxID', 'ScientificName', 'SampleName', 'Sex', 'Tumor', 'Submission',
    'Consent', 'RunHash', 'ReadHash'
])

if args.output:
    output_file_name = args.output
else:
    output_file_name = args.taxid + '.tsv'

with open(output_file_name, 'w') as f:
    f.write(output_header + '\n')
    f.write(proc_xtract.stdout.read().decode('utf-8'))
