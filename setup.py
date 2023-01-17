from urllib.request import urlretrieve
from stat import S_IXUSR, S_IXOTH, S_IXGRP, S_IRUSR, S_IROTH, S_IRGRP, S_IWUSR
from os import mkdir, chmod, remove, chdir
from os.path import dirname, abspath, join, exists, basename
from zipfile import ZipFile
import tarfile
from setuptools import setup, find_packages
import subprocess 

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
    'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip',
    join(lib_dir, 'Trimmomatic-0.38.zip'))

print('Downloading HISAT2 ...')
urlretrieve(
    'https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download',
    join(lib_dir, 'hisat2-2.1.0-Linux_x86_64.zip'))

print('Downloading BBMap ...')
urlretrieve(
    'https://downloads.sourceforge.net/project/bbmap/BBMap_38.00.tar.gz',
    join(lib_dir, 'BBMap_38.00.tar.gz'))

print('Downloading Picard ...')
urlretrieve(
    'https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar',
    join(lib_dir, 'picard.jar'))

print('Downloading GATK v4 ...')
urlretrieve('https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip',
    join(lib_dir, 'gatk-4.1.6.0.zip'))

print('Downloading bam_to_bigwig ...')
subprocess.run(['git', 'clone', 'https://github.com/NAL-i5K/bam_to_bigwig.git', join(lib_dir, 'bam_to_bigwig')])

print('Downloading regtools ...')
subprocess.run(['git', 'clone', 'https://github.com/griffithlab/regtools', join(lib_dir, 'regtools')])

mkdir(join(lib_dir, 'regtools', 'build'))

chdir(join(lib_dir, 'regtools', 'build'))
subprocess.run(['cmake', '..'])
subprocess.run(['make'])
chdir(join(project_root))

print('Unpacking fastQC ...')
with ZipFile(join(lib_dir, 'fastqc_v0.11.7.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)
# fix the permission of fastqc
chmod(
    join(lib_dir, 'FastQC', 'fastqc'),
    S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR | S_IROTH | S_IRGRP | S_IWUSR)

print('Unpacking Trimmomatic ...')
with ZipFile(join(lib_dir, 'Trimmomatic-0.38.zip'), 'r') as zip_ref:
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

print('Unpacking gatk ...')
with ZipFile(join(lib_dir, 'gatk-4.1.6.0.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)

print('Cleaning the files ...')
files = [
    'BBMap_38.00.tar.gz', 'fastqc_v0.11.7.zip',
    'hisat2-2.1.0-Linux_x86_64.zip', 'Trimmomatic-0.38.zip',
    'gatk-4.1.6.0.zip'
]
for f in files:
    remove(join(lib_dir, f))

setup(
    name='rnannot',
    version='0.0.1',
    install_requires=['six'],
    packages=find_packages('.'),
    scripts=['rnannot/RNAseq_annotate.py','rnannot/download_sra_metadata.py'],
    include_package_data=True,
    author='Yi Hsiao',
    author_email='hsiaoyi0504@gmail.com',
    license='Public Domain',
    keywords='RNA-Seq annotation bioinformatics example',
    url=
    'https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline',  # project home page, if any
    )
    
