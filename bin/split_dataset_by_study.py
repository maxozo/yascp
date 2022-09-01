#!/usr/bin/env python3

## use the donor assgnment file produced by gtcheck_assign_summary.py
## to split the dataset for handover by study

import sys
import os
import csv

FOFN_NAME_INPUT = "files.tsv"
OUTDIR_UNKNOWN = "GT_UNRESOLVED"
STUDY_UNKNOWN = "NONE"
NONE_SYMBOL = 'N/A'

FNAM_MANIFEST = "manifest.tsv"
MANIFEST_HEADER = \
    "pool_id\tdonor_id\tstudy\t"\
    "filename_h5ad_encrypted\tchecksum_h5ad_encrypted\t" \
    "filename_annotation_tsv_encrypted\tchecksum_annotation_tsv_encrypted\n"
MANIFEST_FMT = \
    "{0[0]:s}\t{0[1]:s}\t{0[2]:s}\t" \
    "{0[3]:s}\t{0[4]:s}\t" \
    "{0[5]:s}\t{0[6]:s}\n"

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

def read_md5_checksum(fnam):
    chksm = []
    with open(fnam+'.md5', 'r') as infh:
        for line in infh:
            s = line.strip()
            if len(s) > 0:
                chksm.append(s)
    return chksm

def write_manifest(outdir, tranche_name, manifest_arr):
    manifest_dict = {}
    fnam_manifest = tranche_name + "_" + FNAM_MANIFEST
    # write overall manifest and make dict by study
    oufnam = os.path.join(outdir, fnam_manifest)
    oufh = open(oufnam, 'w')
    oufh.write(MANIFEST_HEADER)
    ctr = 0
    for t in manifest_arr:
        ctr += 1
        study_dir = t[2]
        study_symbol = convert_study_dir_to_symbol(study_dir)
        oufh.write(MANIFEST_FMT.format(t[:2]+(study_symbol,)+t[3:]))
        try:
            manifest_dict[study_dir].append(t)
        except KeyError:
            manifest_dict[study_dir] = [t]
    oufh.close()
    # now write per-study manifests
    for study_dir in manifest_dict:
        dirpath = os.path.join(outdir, study_dir)
        if os.path.exists(dirpath):
            study_symbol = convert_study_dir_to_symbol(study_dir)
            fpath_manifest = os.path.join(dirpath, fnam_manifest)
            oufh = open(fpath_manifest, 'w')
            oufh.write(MANIFEST_HEADER)
            for t in manifest_dict[study_dir]:
                oufh.write(MANIFEST_FMT.format(t[:2]+(study_symbol,)+t[3:]))
            oufh.close()
        else:
            sys.stderr.write("WARNING: directory {:s} does not exist.".format(dirpath))
    oufh.close()
    return ctr

def make_symlinks(dst_dir, fnams, src_dir):
    for fn in fnams:
        if fn is None:
            sys.stderr.write("WARNING: missing files in {:s}\n".format(src_dir))
        else:
            src_path = os.path.join(src_dir, fn)
            dst_path = os.path.join(dst_dir, fn)
            print("setting link {:s} -> {:s}".format(src_path, dst_path))
            os.symlink(src_path, dst_path)
    return

def get_file_list_dict(pool_direntry):
    file_lst = {} # key: donor_id, value: [h5ad file, tsv file]
    for fnam in os.listdir(pool_direntry.path):
        fnfld = fnam.split('.')
        if len(fnfld) < 4 or fnfld[-1] != 'gpg':
            continue
        fn_pool_id, donor_id, fn_ext = tuple(fnfld[:3])
        if fn_ext == 'h5ad':
            i = 0
        elif fn_ext == 'tsv':
            i = 1
        elif fn_ext == 'cram':
            i = 3
        else:
            continue
        if donor_id not in file_lst:
            file_lst[donor_id] = [None,None,None]
        file_lst[donor_id][i] = fnam
    return file_lst

def split_dataset(study_dict, dirnam_input, pool_id, dirnam_output):
    # if dirnam_encrypted is not None: fetch engrypted *.gpg and *.gpg.md5 file from that directory
    print("\n=====================\n==+ split_dataset ==+\n=====================")
    manifest_dict = {}
    manifest_arr = []
    # create directories
    for dn in study_dict.keys():
        if dn == STUDY_UNKNOWN:
            d = OUTDIR_UNKNOWN
        else:
            d = dn
        dstn = os.path.join(dirnam_output, d)
        if not os.access(dstn, os.F_OK):
            os.mkdir(dstn)
    #fofn_name = os.path.join(dirnam_input, FOFN_NAME_INPUT)
    #print("opening file ", fofn_name)
    #infh = open(fofn_name, 'r')
    #for row in csv.DictReader(infh, delimiter='\t'):

    for pool_direntry in os.scandir(dirnam_input):
        if not pool_direntry.is_dir() or pool_id != pool_direntry.name:
            continue
        file_lst = get_file_list_dict(pool_direntry)
        for study_name in study_dict:
            if study_name == STUDY_UNKNOWN:
                dstdn = OUTDIR_UNKNOWN
            else:
                dstdn = study_name
            dst_dir = os.path.join(dirnam_output, dstdn)
            for donor_id in study_dict[study_name]:
                src_dir = os.path.join(dirnam_input, pool_id)
                fnams = file_lst[donor_id]
                print(donor_id, fnams)
                make_symlinks(dst_dir, fnams, src_dir)
                t = [pool_id, donor_id, dstdn]
                for fn in fnams:
                    if fn is None:
                        fn = NONE_SYMBOL
                        csm = NONE_SYMBOL
                    else:
                        check_sums = read_md5_checksum(os.path.join(src_dir, fn))
                        if check_sums:
                            csm = check_sums[0]
                        else:
                            csm = NONE_SYMBOL
                            sys.stderr.write("WARNING: could not access {:s}.md5 file.\n".format(fn))
                    t.append(fn)
                    t.append(csm)
                t = tuple(t)
                try:
                    manifest_dict[dstdn].append(t)
                except KeyError:
                    manifest_dict[dstdn] = [t]
                manifest_arr.append(t)
    #infh.close()
    # write_manifests(dirnam_output, pool_id, manifest_dict)
    return manifest_arr

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
    manifest = []
    study_dict, tranche_name = load_donor_assignments_from_assignments_all_pools_tsv(fnam_donor_assignments)
    for pool_id in study_dict:
        manifest_pool = split_dataset(study_dict[pool_id], dirnam_input, pool_id, dirnam_output)
        manifest.extend(manifest_pool)
    write_manifest(dirnam_output, tranche_name, manifest)
    sys.exit(0)
