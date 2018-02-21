from urllib.request import urlretrieve
from os import mkdir
from os.path import dirname, abspath, join
from zipfile import ZipFile

project_root = dirname(abspath(__file__))
software_dir = join(project_root, 'software')

mkdir(software_dir)

print('Downloading fastQC ...')
urlretrieve('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip', join(software_dir, 'fastqc_v0.11.7.zip'))
print('Downloading Trimmomatic ...')
urlretrieve('http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip', join(software_dir, 'Trimmomatic-0.36.zip'))
print('Unpacking fastQC ...')
with ZipFile(join(software_dir, 'fastqc_v0.11.7.zip'), 'r') as zip_ref:
    zip_ref.extractall(join(software_dir,'fastqc_v0.11.7'))
print('Unpacking Trimmomatic ...')
with ZipFile(join(software_dir, 'Trimmomatic-0.36.zip'), 'r') as zip_ref:
    zip_ref.extractall(join(software_dir, 'Trimmomatic-0.36'))

