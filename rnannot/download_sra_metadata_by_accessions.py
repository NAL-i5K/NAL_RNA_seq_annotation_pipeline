import subprocess
import argparse
import time
import re

parser = argparse.ArgumentParser()
parser.add_argument("-a","--accession", help="your accessions",type=str) # replace tax_id
parser.add_argument("-o","--output", help="directory of output folder at, if not specified, use current folder",type=str)                        
args = parser.parse_args()
accession=args.accession  # replace tax_id
print('Processing accession:'+args.accession)  # replace tax_id

proc_esearch = subprocess.Popen(['esearch', '-db', 'sra', '-query' , accession], stdout=subprocess.PIPE)
time.sleep(0.7)
proc_efetch = subprocess.Popen(['efetch', '-format', 'runinfo', '-mode', 'xml'], stdin=proc_esearch.stdout, stdout=subprocess.PIPE)
time.sleep(0.7)
proc_xtract = subprocess.Popen(['xtract', '-pattern', 'Row', '-def', 'N/A','-element', 'Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates', 'avgLength', 'size_MB', 'download_path', 'Experiment', 'LibraryName', 'LibraryStrategy', 'LibrarySelection', 'LibrarySource', 'LibraryLayout', 'InsertSize', 'InsertDev', 'Platform', 'Model', 'SRAStudy', 'BioProject', 'Study_Pubmed_id', 'ProjectID', 'Sample', 'BioSample', 'SampleType', 'TaxID', 'ScientificName', 'SampleName', 'Sex', 'Tumor', 'Submission', 'Consent', 'RunHash', 'ReadHash' ], stdin=proc_efetch.stdout, stdout=subprocess.PIPE)
if args.output:
    with open(args.output, 'w') as f:
         f.write('\t'.join('Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Sex,Tumor,Submission,Consent,RunHash,ReadHash'.split(',')) + '\n')
         f.write(proc_xtract.stdout.read().decode('utf-8') )         
else:
    with open(args.accession+'.tsv', 'w') as f:  # replace tax_id
         f.write('\t'.join('Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Sex,Tumor,Submission,Consent,RunHash,ReadHash'.split(',')) + '\n')
         f.write(proc_xtract.stdout.read().decode('utf-8') )
