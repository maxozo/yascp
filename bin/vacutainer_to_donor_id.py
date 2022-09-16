#!/usr/bin/env python3

## translate donor ids from a conversion table

import sys
import csv

nargs = len(sys.argv)

def load_bridge_file (fnam):
    id_dict = {}
    infh = open(fnam_conversion_table, 'r')
    reader = csv.DictReader(infh, delimiter='\t')
    for r in reader:
        # print(r)
        v_id = r['vacutainer_id']
        d_id = r['vcf_donor_id']
        if d_id.startswith("Unavailable"):
            continue
        if v_id in id_dict and d_id != id_dict[v_id]:
            sys.exit("ERROR: multiple entries with the same key: {:s}: {:s}, {:s}\n".format(v_id, id_dict[v_id], d_id))
        id_dict[v_id] = d_id
    infh.close()
    return id_dict

if __name__ == '__main__':
    if nargs != 4:
        sys.exit("usage: {:s} <conversion table> <comma separated list of ids> <output file>"
            .format(sys.argv[0]))

    fnam_conversion_table = sys.argv[1]
    id_list_str = sys.argv[2]
    fnam_output = sys.argv[3]

    id_dict = load_bridge_file(fnam_conversion_table)

    oufh = open(fnam_output, 'w')
    vacutainer_ids = id_list_str.split(",")
    print("vacutainer_ids", vacutainer_ids)

    for vacu_id in vacutainer_ids:
        if vacu_id not in id_dict:
            vacu_id = vacu_id.lstrip('0') # bridge file/VCF don't contain leading 0s
        try:
            donor_id = id_dict[vacu_id]
        except KeyError:
            pass
        else:
            oufh.write(donor_id + '\n')

    oufh.close()
    sys.exit(0)
