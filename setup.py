from urllib.request import urlretrieve
from stat import S_IXUSR, S_IXOTH, S_IXGRP, S_IRUSR, S_IROTH, S_IRGRP, S_IWUSR
from os import mkdir, chmod
from os.path import dirname, abspath, join
from zipfile import ZipFile
import tarfile


project_root = dirname(abspath(__file__))
lib_dir = join(project_root, 'rnannot', 'lib')

mkdir(lib_dir)

print('Downloading fastQC ...')
urlretrieve('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip', join(lib_dir, 'fastqc_v0.11.7.zip'))

print('Downloading Trimmomatic ...')
urlretrieve('http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip', join(lib_dir, 'Trimmomatic-0.36.zip'))

print('Downloading HISAT2 ...')
urlretrieve('ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip', join(lib_dir, 'hisat2-2.1.0-Linux_x86_64.zip'))

print('Downloading BBMap ...')
urlretrieve('https://downloads.sourceforge.net/project/bbmap/BBMap_38.00.tar.gz', join(lib_dir, 'BBMap_38.00.tar.gz'))

print('Unpacking fastQC ...')
with ZipFile(join(lib_dir, 'fastqc_v0.11.7.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)
# fix the permission of fastqc
chmod(join(lib_dir, 'FastQC', 'fastqc'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)

print('Unpacking Trimmomatic ...')
with ZipFile(join(lib_dir, 'Trimmomatic-0.36.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)

print('Unpacking HISAT2 ...')
with ZipFile(join(lib_dir, 'hisat2-2.1.0-Linux_x86_64.zip'), 'r') as zip_ref:
    zip_ref.extractall(lib_dir)
# fix the permission of HISAT2
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2-align-s'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2-align-l'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2-build'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2-build-s'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)
chmod(join(lib_dir, 'hisat2-2.1.0', 'hisat2-build-l'), S_IXUSR | S_IXOTH | S_IXGRP | S_IRUSR  | S_IROTH | S_IRGRP | S_IWUSR)

print('Unpacking BBMap ...')
tar = tarfile.open(join(lib_dir, 'BBMap_38.00.tar.gz'), 'r:gz')
tar.extractall(lib_dir)
tar.close()

# TODO: clean up files

