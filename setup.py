from setuptools import setup, find_packages

setup(
    name='rnannot',
    version='2026.10.26',
    python_requires='>=3.10',
    install_requires=[],          # add 'six' back only if the scripts import it
    packages=find_packages('.'),
    scripts=['rnannot/RNAseq_annotate.py', 'rnannot/download_sra_metadata.py'],
    include_package_data=True,
    author='Yi Hsiao',
    author_email='hsiaoyi0504@gmail.com',
    license='Public Domain',
    keywords='RNA-Seq annotation bioinformatics',
    url='https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline',
)
