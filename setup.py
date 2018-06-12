from urllib.request import urlretrieve
from stat import S_IXUSR, S_IXOTH, S_IXGRP, S_IRUSR, S_IROTH, S_IRGRP, S_IWUSR
from os import mkdir, chmod, remove
from os.path import dirname, abspath, join, exists, basename
from zipfile import ZipFile
import tarfile
from setuptools import setup, find_packages

project_root = dirname(abspath(__file__))
lib_dir = join(project_root, 'rnannot', 'lib')

if not exists(lib_dir):
    mkdir(lib_dir)

print('Downloading fastQC ...')
urlretrieve(
    'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip',
    join(lib_dir, 'fastqc_v0.11.7.zip'))

print('Downloading Trimmomatic ...')
urlretrieve(
    'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip',
    join(lib_dir, 'Trimmomatic-0.36.zip'))

print('Downloading HISAT2 ...')
urlretrieve(
    'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip',
    join(lib_dir, 'hisat2-2.1.0-Linux_x86_64.zip'))

print('Downloading BBMap ...')
urlretrieve(
    'https://downloads.sourceforge.net/project/bbmap/BBMap_38.00.tar.gz',
    join(lib_dir, 'BBMap_38.00.tar.gz'))

print('Downloading Picard ...')
urlretrieve(
    'https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar',
    join(lib_dir, 'picard.jar'))

print('Downloading GATK v3 ...')
urlretrieve('https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef',
    join(lib_dir, 'GenomeAnalysisTK.tar.bz2'))

print('Unpacking fastQC ...')
with ZipFile(join(lib_dir, 'fastqc_v0.11.7.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)
# fix the permission of fastqc
chmod(
    join(lib_dir, 'FastQC', 'fastqc'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)

print('Unpacking Trimmomatic ...')
with ZipFile(join(lib_dir, 'Trimmomatic-0.36.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)

print('Unpacking HISAT2 ...')
with ZipFile(join(lib_dir, 'hisat2-2.1.0-Linux_x86_64.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)
# fix the permission of HISAT2
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2-align-s'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2-align-l'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2-build'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2-build-s'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)
chmod(
    join(lib_dir, 'hisat2-2.1.0', 'hisat2-build-l'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)

print('Unpacking BBMap ...')
tar = tarfile.open(join(lib_dir, 'BBMap_38.00.tar.gz'), 'r:gz')
tar.extractall(lib_dir)
tar.close()

print('Unpacking GATK v3 ...')
tar = tarfile.open(join(lib_dir, 'GenomeAnalysisTK.tar.bz2'), 'r:bz2')
for member in tar.getmembers():
    if member.isreg():  # skip if the TarInfo is not files
        member.name = basename(member.name) # remove the path by reset it
        tar.extract(member, lib_dir) # extract
tar.close()

print('Cleaning the files ...')
files = [
    'BBMap_38.00.tar.gz', 'fastqc_v0.11.7.zip',
    'hisat2-2.1.0-Linux_x86_64.zip', 'Trimmomatic-0.36.zip',
    'GenomeAnalysisTK.tar.bz2'
]
for f in files:
    remove(join(lib_dir, f))

setup(
    name='rnannot',
    version='0.0.1',
    packages=find_packages(),
    scripts=['rnannot/RNAseq_annotate.py'],
    include_package_data=True,
    author='Yi Hsiao',
    author_email='hsiaoyi0504@gmail.com',
    license='Public Domain',
    keywords='RNA-Seq annotation bioinformatics example',
    url=
    'https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline',  # project home page, if any
    )
