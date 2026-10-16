from os import environ, path

# Tool versions. Defaults match the Dockerfile ARGs; override with
# environment variables so future upgrades only need a Dockerfile change.
TRIMMOMATIC_VERSION = environ.get('TRIMMOMATIC_VERSION', '0.39')
HISAT2_VERSION = environ.get('HISAT2_VERSION', '2.2.1')


def get_lib_path():
    # RNANNOT_LIB lets the tools live outside the package directory if needed.
    return environ.get(
        'RNANNOT_LIB',
        path.join(path.abspath(path.dirname(__file__)), 'lib'))


def get_trimmomatic_jar_path():
    return path.join(get_lib_path(), 'Trimmomatic-' + TRIMMOMATIC_VERSION,
                     'trimmomatic-' + TRIMMOMATIC_VERSION + '.jar')


def get_fastqc_path():
    return path.join(get_lib_path(), 'FastQC', 'fastqc')


def get_trimmomatic_adapter_path(file_name):
    return path.join(get_lib_path(), 'Trimmomatic-' + TRIMMOMATIC_VERSION,
                     'adapters', file_name)


def get_bbmap_adapter_path():
    return path.join(get_lib_path(), 'bbmap', 'resources', 'adapters.fa')


def get_hisat2_command_path(cmd):
    return path.join(get_lib_path(), 'hisat2-' + HISAT2_VERSION, cmd)


def get_bbmap_command_path(cmd):
    return path.join(get_lib_path(), 'bbmap', cmd)


def get_picard_jar_path():
    return path.join(get_lib_path(), 'picard.jar')


def get_bam_to_bigwig_path():
    # Unmodified script from the NAL-i5K/bam_to_bigwig repo, cloned into lib/
    return path.join(get_lib_path(), 'bam_to_bigwig', 'bam_to_bigwig.py')


def get_regtools_path():
    return path.join(get_lib_path(), 'regtools', 'build', 'regtools')
