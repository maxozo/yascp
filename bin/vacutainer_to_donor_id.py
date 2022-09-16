#!/usr/bin/env python3

## translate donor ids from a conversion table

import sys
nargs = len(sys.argv)

if nargs != 4:
    sys.exit("usage: {:s} <conversion table> <comma separated list of ids> <output file>"
        .format(sys.argv[0]))

fnam_conversion_table = sys.argv[1]
id_list_str = sys.argv[2]
fnam_output = sys.argv[3]

oufh = open(fnam_output, 'w')
vacutainer_ids = id_list_str.split(",")
print("vacutainer_ids", vacutainer_ids)
with open(fnam_conversion_table, 'r') as infh:
    for lin in infh:
        fld = lin.split()
        print(fld)

        if fld[1] in vacutainer_ids or \
            fld[1].lstrip('0') in vacutainer_ids:
            # bridge file doesn't contain leading 0s for UKBB IDs like '0030007452847'
            oufh.write("{:s}\n".format(fld[0]))

oufh.close()
sys.exit(0)
