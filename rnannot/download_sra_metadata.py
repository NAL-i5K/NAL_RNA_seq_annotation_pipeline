import subprocess
import argparse
import time
import re

parser = argparse.ArgumentParser()
parser.add_argument("-t","--taxid", help="your tax_id",type=str)
parser.add_argument("-o","--output", help="directory of output folder at, if not specified, use current folder",type=str)                        
args = parser.parse_args()
tax_id=args.taxid
print('Processing tax id:'+args.taxid)

proc_esearch = subprocess.Popen(['esearch', '-db', 'taxonomy', '-query' , tax_id + '[uid]'], stdout=subprocess.PIPE)
time.sleep(0.7)
proc_elink = subprocess.Popen(['elink', '-target', 'sra'], stdin=proc_esearch.stdout, stdout=subprocess.PIPE)
time.sleep(0.7)
proc_efilter = subprocess.Popen(['efilter', '-query', 'rna seq[stra] AND Transcriptomic[src] NOT PACBIO_SMRT[platform]'], stdin=proc_elink.stdout, stdout=subprocess.PIPE)
time.sleep(0.7)
proc_efetch = subprocess.Popen(['efetch', '-format', 'runinfo', '-mode', 'xml'], stdin=proc_efilter.stdout, stdout=subprocess.PIPE)
time.sleep(0.7)
proc_xtract = subprocess.Popen(['xtract', '-pattern', 'Row', '-def', 'N/A','-element', 'Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates', 'avgLength', 'size_MB', 'download_path', 'Experiment', 'LibraryName', 'LibraryStrategy', 'LibrarySelection', 'LibrarySource', 'LibraryLayout', 'InsertSize', 'InsertDev', 'Platform', 'Model', 'SRAStudy', 'BioProject', 'Study_Pubmed_id', 'ProjectID', 'Sample', 'BioSample', 'SampleType', 'TaxID', 'ScientificName', 'SampleName', 'Sex', 'Tumor', 'Submission', 'Consent', 'RunHash', 'ReadHash' ], stdin=proc_efetch.stdout, stdout=subprocess.PIPE)
if args.output:
    with open(args.output, 'w') as f:
         f.write('\t'.join('Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Sex,Tumor,Submission,Consent,RunHash,ReadHash'.split(',')) + '\n')
         f.write(proc_xtract.stdout.read().decode('utf-8') )         
else:
    with open(args.taxid+'.tsv', 'w') as f:
         f.write('\t'.join('Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Sex,Tumor,Submission,Consent,RunHash,ReadHash'.split(',')) + '\n')
         f.write(proc_xtract.stdout.read().decode('utf-8') )
