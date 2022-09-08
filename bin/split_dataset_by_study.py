#!/usr/bin/env python3

## use the donor assgnment file produced by gtcheck_assign_summary.py
## to split the dataset for handover by study

import sys
import os
import re
import csv
import copy

FOFN_NAME_INPUT = "files.tsv"
OUTDIR_UNKNOWN = "GT_UNRESOLVED"
STUDY_UNKNOWN = "NONE"
NONE_SYMBOL = 'N/A'

FNAM_MANIFEST = "manifest.tsv"
MANIFEST_HEADER = \
    "pool_id\tdonor_id\tstudy\t"\
    "filename_h5ad_encrypted\tchecksum_h5ad_encrypted\tchecksum_h5ad\t" \
    "filename_annotation_tsv_encrypted\tchecksum_annotation_tsv_encrypted\tchecksum_annotation_tsv\t" \
    "filename_cram_encrypted\tchecksum_cram_encrypted\tchecksum_cram\n"
MANIFEST_FMT = \
    "{0[0]:s}\t{0[1]:s}\t{0[2]:s}\t" \
    "{1[0]:s}\t{1[1]:s}\t{1[2]:s}\t" \
    "{2[0]:s}\t{2[1]:s}\t{2[2]:s}\t" \
    "{3[0]:s}\t{3[1]:s}\t{3[2]:s}\n"

FNAM_POOL_DONOR_MATCH = re.compile("^(\S+)_barcodes_(\S+)_possorted_bam$")

FILE_EXTENSIONS = ('.h5ad', '.tsv', '.cram') # don't change the order
FILE_EXT_ENCRYPTED = '.gpg'
FILE_EXT_CHECKSUM = '.md5'

def convert_study_dir_to_symbol(study_dir):
    if study_dir == OUTDIR_UNKNOWN:
        study_symbol = NONE_SYMBOL
    elif study_dir.startswith('GT_'):
        study_symbol = study_dir[3:]
    else:
        study_symbol = study_dir
    return study_symbol

def load_donor_assignments(infn):
    study_dict = {}
    with open(infn, 'r') as infh:
        for row in csv.DictReader(infh):
            study_name = row['panel']
            donor_id = row['donor_query']
            try:
                study_dict[study_name].append(donor_id)
            except KeyError:
                study_dict[study_name] = [donor_id]
    return study_dict

def load_donor_assignments_from_assignments_all_pools_tsv(fnam):
    pools_dict = {}
    tranche_name = None
    with open(fnam, 'r') as infh:
        for row in csv.DictReader(infh, delimiter='\t'):
            donor_id = row['donor_query']
            donor_panel = row['final_panel']
            pool_id = row['pool']
            if tranche_name is None:
                tranche_name = row['tranche']
            elif tranche_name != row['tranche']:
                sys.exit(
                    "ERROR: encountered multiple tranche names: {:s}, {:s}\n"
                    .format(tranche_name, row['tranche'])
                )
            if pool_id not in pools_dict:
                pools_dict[pool_id] = {donor_panel: [donor_id]}
            elif donor_panel not in pools_dict[pool_id]:
                pools_dict[pool_id][donor_panel] = [donor_id]
            else:
                pools_dict[pool_id][donor_panel].append(donor_id)
    return pools_dict, tranche_name

def load_dirs_from_file(fnam):
    dirnams = []
    with open(fnam, 'r') as infh:
        for line in infh:
            dns = line.split()
            dirnams.extend(dns)
    return dirnams

def read_md5_checksum(fnam):
    chksm = []
    with open(fnam+FILE_EXT_CHECKSUM, 'r') as infh:
        for line in infh:
            s = line.strip()
            if len(s) > 0:
                chksm.append(s)
    return chksm

def read_md5_checksum_unencrypted(fnam):
    filnam_ext = os.path.splitext(fnam)
    if filnam_ext[1] != FILE_EXT_ENCRYPTED:
        sys.exit("ERROR: unexpected filename '{:s}' doesn't have extension {:s}"
            .format(fnam, FILE_EXT_ENCRYPTED))
    return read_md5_checksum(filnam_ext[0])

def merge_manifest_rows(manifest_pool, manifest_crams):
    manifest = []
    for donor_id in manifest_pool:
        manifest_row = copy.deepcopy(manifest_pool[donor_id])
        if manifest_crams:
            row2 = manifest_crams[donor_id]
            if manifest_row[:3] != row2[:3]:
                sys.exit("ERROR: manifest mismatch.")
            manifest_row[5] = row2[5]
        manifest.append(manifest_row)
    return manifest

def write_manifest(outdir, tranche_name, manifest_arr, fnam_manifest = FNAM_MANIFEST):
    manifest_dict = {}
    fnam_manifest = tranche_name + "_" + fnam_manifest
    # write overall manifest and make dict by study
    oufnam = os.path.join(outdir, fnam_manifest)
    oufh = open(oufnam, 'w')
    oufh.write(MANIFEST_HEADER)
    ctr = 0
    for row in manifest_arr:
        ctr += 1
        study_dir = row[2]
        row[2] = convert_study_dir_to_symbol(study_dir)
        print("data_tup", row)
        oufh.write(MANIFEST_FMT.format(row[:3], row[3], row[4], row[5]))
        try:
            manifest_dict[study_dir].append(row)
        except KeyError:
            manifest_dict[study_dir] = [row]
    oufh.close()
    # now write per-study manifests
    for study_dir in manifest_dict:
        dirpath = os.path.join(outdir, study_dir)
        if os.path.exists(dirpath):
            fpath_manifest = os.path.join(dirpath, fnam_manifest)
            oufh = open(fpath_manifest, 'w')
            oufh.write(MANIFEST_HEADER)
            for row in manifest_dict[study_dir]:
                oufh.write(MANIFEST_FMT.format(row[:3], row[3], row[4], row[5]))
            oufh.close()
        else:
            sys.stderr.write("WARNING: directory {:s} does not exist.".format(dirpath))
    oufh.close()
    return ctr

def make_symlinks(dst_dir, fpaths):
    for fpath in fpaths:
        #if fn is None:
        #    sys.stderr.write("WARNING: missing files in {:s}\n".format(src_dir))
        #else:
        if fpath is not None:
            dst_path = os.path.join(dst_dir, os.path.basename(fpath))
            print("setting link {:s} -> {:s}".format(fpath, dst_path))
            os.symlink(fpath, dst_path)
    return

def get_file_list_dict(pool_direntry):
    file_lst = {} # key: donor_id, value: {file_extension: file_name}]
    for fnam in os.listdir(pool_direntry.path):
        donor_id = None
        i = -1
        fn, ext = os.path.splitext(fnam)
        print("os.patwh.splitext: {0:s}, {1:s}, {2:s}".format(fnam, fn, ext))
        if ext == FILE_EXT_ENCRYPTED:
            fn_base, fn_typ = os.path.splitext(fn)
            try:
                i = FILE_EXTENSIONS.index(fn_typ)
            except ValueError:
                pass
            else:
                fnfld = fn_base.split('.')
                print("get_file_list_dict", fnam, fnfld)
                if len(fnfld) == 2:
                    fn_pool_id, donor_id = tuple(fnfld)
                else:
                    mm = FNAM_POOL_DONOR_MATCH.match(fn_base)
                    if mm:
                        fn_pool_id, donor_id = mm.groups()
                if i >= 0 and donor_id is not None:
                    if donor_id not in file_lst:
                        file_lst[donor_id] = [None,None,None]
                    file_lst[donor_id][i] = os.path.join(pool_direntry.path, fnam)
    return file_lst

def get_file_paths_dict(dirnam_input, pool_id):
    filpaths = {}
    is_list = type(dirnam_input) is list
    if is_list:
        input_dir = os.curdir
    else:
        input_dir= dirnam_input
    for direntry in os.scandir(input_dir):
        if direntry.is_dir() and direntry.name.startswith(pool_id) \
            and ((not is_list) or direntry.name in dirnam_input):
            filpaths.update(get_file_list_dict(direntry))
    return filpaths

def get_manifest_entries_for_encrytped_file_path(fpath):
    csm = NONE_SYMBOL
    csm_unencrypted = NONE_SYMBOL
    if fpath is None:
        fn = NONE_SYMBOL
    else:
        fn = os.path.basename(fpath)
        check_sums = read_md5_checksum(fpath)
        if check_sums:
            csm = check_sums[0]
        else:
            sys.stderr.write("WARNING: could not access {:s}{:s} file.\n".format(fn,FILE_EXT_CHECKSUM))
        check_sums_unencrypted = read_md5_checksum_unencrypted(fpath)
        if check_sums_unencrypted:
            csm_unencrypted = check_sums_unencrypted[0]
        else:
            sys.stderr.write("WARNING: could not access unencrypted {:s} checksum file.\n".format(fn))
    return fn, csm, csm_unencrypted

def split_dataset(study_dict, dirnam_input, pool_id, dirnam_output):
    # if dirnam_encrypted is not None: fetch engrypted *.gpg and *.gpg.md5 file from that directory
    print("\n=====================\n==+ split_dataset ==+\n=====================")
    manifest_dict = {}
    # create directories
    for dn in study_dict.keys():
        if dn == STUDY_UNKNOWN:
            d = OUTDIR_UNKNOWN
        else:
            d = dn
        dstn = os.path.join(dirnam_output, d)
        if not os.access(dstn, os.F_OK):
            os.mkdir(dstn)
    filpath_dict = get_file_paths_dict(dirnam_input, pool_id)
    print("filpath_dict", filpath_dict)
    for study_name in study_dict:
        if study_name == STUDY_UNKNOWN:
            dstdn = OUTDIR_UNKNOWN
        else:
            dstdn = study_name
        dst_dir = os.path.join(dirnam_output, dstdn)
        for donor_id in study_dict[study_name]:
            fpaths = filpath_dict[donor_id]
            print(donor_id, fpaths)
            make_symlinks(dst_dir, fpaths)
            manifest_row = [pool_id, donor_id, dstdn]
            for fpath in fpaths:
                data_tup = get_manifest_entries_for_encrytped_file_path(fpath)
                manifest_row.append(data_tup)
            if donor_id in manifest_dict:
                sys.exit("ERROR: donor '{:s}' occurs multiple time in pool '{:s}'"
                    .format(donor_id, pool_id))
            manifest_dict[donor_id] = manifest_row
    return manifest_dict

if __name__ == '__main__':
    nargs = len(sys.argv)

    if nargs < 4 or nargs > 5:
        sys.exit(
            "usage: {:s} <donor assignment file [tsv]> <input directory> <output_directory> [<CRAM dirs list>]"
            .format(sys.argv[0])
        )

    fnam_donor_assignments = sys.argv[1]
    dirnam_input = sys.argv[2]
    dirnam_output = sys.argv[3]
    cram_dirs_input = None
    if nargs > 4:
        cram_dirs_input = sys.argv[4]
    os.makedirs(dirnam_output, exist_ok = False)

    #study_dict = load_donor_assignments(fnam_donor_assignments)
    manifest_rows = []
    cram_dirs = []
    if cram_dirs_input:
        cram_dirs = load_dirs_from_file(cram_dirs_input)
    print("CRAM directories:",cram_dirs)
    study_dict, tranche_name = load_donor_assignments_from_assignments_all_pools_tsv(fnam_donor_assignments)
    for pool_id in study_dict:
        manifest_pool = split_dataset(study_dict[pool_id], dirnam_input, pool_id, dirnam_output)
        print("manifest_pool", manifest_pool)
        if cram_dirs:
            manifest_crams = split_dataset(study_dict[pool_id], cram_dirs, pool_id, dirnam_output)
            print("manifest_crams", manifest_crams)
        else:
            manifest_crams = None
        manifest_joined_rows = merge_manifest_rows(manifest_pool, manifest_crams)
        print("manifest_joined_rows", manifest_joined_rows)
        manifest_rows.extend(manifest_joined_rows)

    write_manifest(dirnam_output, tranche_name, manifest_rows, fnam_manifest = FNAM_MANIFEST)
    sys.exit(0)
