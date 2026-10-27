"""Shared i5k naming rules for add_trackList.py and move_data.py."""
import re
from os import path

# Rank words skipped when parsing scientific names
RANK_MARKERS = {'subsp.', 'subsp', 'ssp.', 'ssp', 'var.', 'var', 'subvar.',
                'f.', 'forma', 'cf.', 'aff.', 'x', '×'}


def species_dir_type(value):
    """argparse type: letters, digits, underscores; case kept."""
    import argparse
    value = value.strip()
    if not re.fullmatch(r'[A-Za-z0-9_]+', value):
        raise argparse.ArgumentTypeError(
            'must contain only letters, digits and underscores, got "{}"'.format(value))
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


def full_name(parts):
    """New convention: Genus_species[_subspecies]. Also the pipeline's file prefix."""
    return '_'.join(parts)


def legacy_codes(parts):
    """Every legacy code a species directory might have been created with."""
    codes = [(parts[0][:3] + parts[1][:3]).lower()]            # original gggsss rule
    if len(parts) == 3:
        codes.append(''.join(p[:2] for p in parts).lower())     # ggssii trinomial rule
    return list(dict.fromkeys(codes))                           # dedupe, keep order


def resolve_species_dir(base_path, scientific_name, override=None):
    """Return (species_dir, reason). Raises SystemExit if the choice is ambiguous."""
    if override:
        return override, 'supplied with --species_dir'
    parts = parse_scientific_name(scientific_name)
    new = full_name(parts)
    legacy_found = [c for c in legacy_codes(parts) if path.isdir(path.join(base_path, c))]
    new_found = path.isdir(path.join(base_path, new))

    if new_found and legacy_found:
        raise SystemExit(
            'Both legacy ({}) and new ({}) directories exist under {}.\n'
            'Choose one with --species_dir.'.format(', '.join(legacy_found), new, base_path))
    if len(legacy_found) > 1:
        raise SystemExit(
            'More than one legacy directory matches: {}. Choose one with --species_dir.'.format(
                ', '.join(legacy_found)))
    if legacy_found:
        return legacy_found[0], 'existing legacy gggsss directory'
    if new_found:
        return new, 'existing full-name directory'
    return new, 'new species: full-name convention'


def analyses_folder(species_dir):
    """Undated analyses folder, shared by all RNA-Seq datasets for an assembly."""
    return species_dir + '-RNA-Seq_v1.0'
