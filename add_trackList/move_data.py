import subprocess
import argparse
import re
import shlex
from os import path

# Rank words skipped when parsing scientific names (same rules as the pipeline and loader)
RANK_MARKERS = {'subsp.', 'subsp', 'ssp.', 'ssp', 'var.', 'var', 'subvar.',
                'f.', 'forma', 'cf.', 'aff.', 'x', '×'}


def gggsss_type(value):
    """Validate a user-supplied organism code: letters/digits only, stored lowercase."""
    value = value.strip().lower()
    if not re.fullmatch(r'[a-z0-9]+', value):
        raise argparse.ArgumentTypeError(
            'must contain only letters and digits, got "{}"'.format(value))
    return value


def parse_scientific_name(scientific_name):
    """Return [genus, species] or [genus, species, subspecies], cleaned for file names."""
    parts = [p for p in scientific_name.split() if p.lower() not in RANK_MARKERS]
    parts = [re.sub(r'[^A-Za-z0-9]', '', p) for p in parts]
    parts = [p for p in parts if p]
    if len(parts) < 2:
        raise ValueError('Cannot parse genus and species from "{}"'.format(scientific_name))
    if len(parts) > 3:
        parts = parts[:2] + [parts[-1]]   # keep the lowest infraspecific epithet
    return parts


def make_gggsss(name_parts):
    """gggsss for binomials, ggssii (2 letters each) for trinomials; lowercase."""
    if len(name_parts) == 3:
        return ''.join(p[:2] for p in name_parts).lower()
    return (name_parts[0][:3] + name_parts[1][:3]).lower()


parser = argparse.ArgumentParser()
parser.add_argument("-Node1a", "--node1_account", required=True,
                    help="apollo-node1 account e.g user@apollo-node1.nal.usda.gov", type=str)
parser.add_argument("-s", "--Source", required=True, help="Source.txt file path", type=str)
parser.add_argument("-g", "--gggsss", type=gggsss_type, default=None,
                    help="organism code; use the same value given to the pipeline and loader "
                         "(default: derived from the scientific name)")
parser.add_argument("--label", default=None,
                    help="track label used by the loader (the BAM file name without .bam). "
                         "Only needed if the pipeline was run with -R/--Run_prefix.")
parser.add_argument("--base_path", default="/app/data/other_species",
                    help="species data root on both servers. Default: %(default)s")
args = parser.parse_args()


def run(cmd):
    """Run a local command and stop if it fails."""
    print('RUN: ' + ' '.join(cmd))
    subprocess.run(cmd, check=True)


def remote(command):
    """Run a shell command string on apollo-node1 and stop if it fails."""
    print('REMOTE: ' + command)
    subprocess.run(['ssh', args.node1_account, command], check=True)


q = shlex.quote

#get info from Source.txt
with open(args.Source) as f:
    scientific_name = f.readline().rstrip('\n')
    assembly_name = f.readline().rstrip('\n')
    bam_file_date = f.readline().rstrip('\n')
    Submission = [line.rstrip('\n') for line in f if line.strip()]

#names must match what the loader created
if args.gggsss:
    gggsss = args.gggsss
    print('Using supplied organism code: {}'.format(gggsss))
else:
    gggsss = make_gggsss(parse_scientific_name(scientific_name))
    print('Derived organism code: {}'.format(gggsss))
folder_name = gggsss + '-RNA-Seq_' + bam_file_date + '_v1.0'
label = args.label or (gggsss + '_' + assembly_name + '_RNA-Seq-alignments_' + bam_file_date)
junctions_folder_name = label + '_junctions'

base = path.join(args.base_path, gggsss, assembly_name)
tracklist = path.join(base, 'jbrowse/data/trackList.json')
analyses_dir = path.join(base, 'scaffold', 'analyses')
new_dir_path = path.join(analyses_dir, folder_name)
junctions_dir = path.join(base, 'jbrowse/data/tracks', junctions_folder_name)
link_path = path.join(base, 'jbrowse/data/analyses')

#check that the local files exist before touching the remote server
for p in (tracklist, new_dir_path, junctions_dir):
    if not path.exists(p):
        raise SystemExit('Not found locally: {}\n'
                         'Check --gggsss / --label against the loader run.'.format(p))

#back up trackList on apollo-node1 server (skip if it doesn't exist yet)
remote('test ! -e {0} || cp {0} {1}'.format(q(tracklist), q(tracklist + '.bak')))

#create new directory on apollo-node1 server
remote('mkdir -p {}'.format(q(new_dir_path)))

#rsync analysis files (bam, csi, bigwig, bed) to new directory on apollo-node1 server
run(['rsync', '-P', '-r', new_dir_path + '/', args.node1_account + ':' + new_dir_path + '/'])

#rsync junction track data to apollo-node1 server
remote('mkdir -p {}'.format(q(path.join(base, 'jbrowse/data/tracks'))))
run(['rsync', '-P', '-r', junctions_dir, args.node1_account + ':' + path.join(base, 'jbrowse/data/tracks') + '/'])

#create the jbrowse/data/analyses symlink only if it doesn't already exist
remote('test -e {0} || ln -s {1} {0}'.format(q(link_path), q(analyses_dir)))

#copy trackList.json last, so it never references tracks that haven't arrived
run(['scp', tracklist, args.node1_account + ':' + tracklist])

print('Finished syncing {} to {}'.format(folder_name, args.node1_account))
