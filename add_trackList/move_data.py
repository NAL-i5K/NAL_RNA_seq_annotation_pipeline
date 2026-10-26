import subprocess
import argparse
import shlex
from os import path

from i5k_naming import (species_dir_type, parse_scientific_name, full_name,
                        resolve_species_dir, analyses_folder)

parser = argparse.ArgumentParser()
parser.add_argument("-Node1a", "--node1_account", required=True,
                    help="apollo-node1 account e.g user@apollo-node1.nal.usda.gov", type=str)
parser.add_argument("-s", "--Source", required=True, help="Source.txt file path", type=str)
parser.add_argument("--species_dir", "-g", "--gggsss", dest="species_dir", type=species_dir_type,
                    default=None,
                    help="species directory (default: same detection as add_trackList.py)")
parser.add_argument("--prefix", type=species_dir_type, default=None,
                    help="file name prefix used by the pipeline (default: Genus_species[_subspecies])")
parser.add_argument("--label", default=None,
                    help="track label used by the loader (the BAM file name without .bam). "
                         "Only needed if the pipeline was run with -R/--Run_prefix.")
parser.add_argument("--base_path", default="/app/data/other_species",
                    help="species data root on both servers. Default: %(default)s")
parser.add_argument("--dry_run", action="store_true",
                    help="print the resolved paths and exit without contacting apollo-node1")
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

#resolve against the local tree; the remote tree mirrors it
species_dir, reason = resolve_species_dir(args.base_path, scientific_name, args.species_dir)
folder_name = analyses_folder(species_dir)
prefix = args.prefix or full_name(parse_scientific_name(scientific_name))
label = args.label or (prefix + '_' + assembly_name + '_RNA-Seq-alignments_' + bam_file_date)
junctions_folder_name = label + '_junctions'

base = path.join(args.base_path, species_dir, assembly_name)
tracklist = path.join(base, 'jbrowse/data/trackList.json')
analyses_dir = path.join(base, 'scaffold', 'analyses')
new_dir_path = path.join(analyses_dir, folder_name)
junctions_dir = path.join(base, 'jbrowse/data/tracks', junctions_folder_name)
link_path = path.join(base, 'jbrowse/data/analyses')

print('Species directory: {} ({})'.format(species_dir, reason))
print('Analyses folder:   {}'.format(new_dir_path))
print('Junctions track:   {}'.format(junctions_dir))

#check that the local files exist before touching the remote server
for p in (tracklist, new_dir_path, junctions_dir):
    if not path.exists(p):
        raise SystemExit('Not found locally: {}\n'
                         'Check --species_dir / --prefix / --label against the loader run.'.format(p))

if args.dry_run:
    print('Dry run: nothing synced.')
    raise SystemExit(0)

#back up trackList on apollo-node1 server (skip if it doesn't exist yet)
remote('test ! -e {0} || cp {0} {1}'.format(q(tracklist), q(tracklist + '.bak')))

#create the analyses folder on apollo-node1 server
remote('mkdir -p {}'.format(q(new_dir_path)))

#rsync the analyses folder; files from earlier datasets are already there and are skipped
run(['rsync', '-P', '-r', new_dir_path + '/', args.node1_account + ':' + new_dir_path + '/'])

#rsync junction track data to apollo-node1 server
remote('mkdir -p {}'.format(q(path.join(base, 'jbrowse/data/tracks'))))
run(['rsync', '-P', '-r', junctions_dir, args.node1_account + ':' + path.join(base, 'jbrowse/data/tracks') + '/'])

#create the jbrowse/data/analyses symlink only if it doesn't already exist
remote('test -e {0} || ln -s {1} {0}'.format(q(link_path), q(analyses_dir)))

#copy trackList.json last, so it never references tracks that haven't arrived
run(['scp', tracklist, args.node1_account + ':' + tracklist])

print('Finished syncing {} to {}'.format(folder_name, args.node1_account))
