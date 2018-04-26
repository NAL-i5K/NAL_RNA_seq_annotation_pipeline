from os import path


def get_trimmomatic_jar_path():
    current_dir = path.abspath(path.dirname(__file__))
    return path.join(current_dir, 'lib', 'Trimmomatic-0.36', 'trimmomatic-0.36.jar')


def get_fastqc_path():
    current_dir = path.abspath(path.dirname(__file__))
    return path.join(current_dir, 'lib', 'FastQC', 'fastqc')


def get_adapter_path(file_name):
    current_dir = path.abspath(path.dirname(__file__))
    return path.join(current_dir, 'lib', 'Trimmomatic-0.36', 'adapters', file_name)

def get_hisat2_command_path(cmd):
    current_dir = path.abspath(path.dirname(__file__))
    return path.join(current_dir, 'lib', 'hisat2-2.1.0', cmd)
