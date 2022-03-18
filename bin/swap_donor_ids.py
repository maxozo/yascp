#!/usr/bin/env python3

## translate donor ids from a conversion table

import sys
import csv

def load_conversion_table(fnam):
    iddict = {}
    with open(fnam, 'r') as infh:
        reader = csv.DictReader(infh, delimiter='\t')
        for row in reader:
            iddict[row['s00046_id']] = row['oragene_id']
    return iddict

def load_full_vcf_ids(fnam):
    donoridd = {}
    linctr = 0
    with open(fnam, 'r') as infh:
        for lin in infh:
            linctr += 1
            fld = lin.split()
            if len(fld) > 1:
                sys.exit("ERROR: expected 1 id per line but got {:d} in file {:s}."
                    .format(len(fld), fnam))
            segm = fld[0].split('_')
            if segm[0] in donoridd:
                if donoridd[segm[0]] != fld[0]:
                    sys.exit("ERROR: donor identifier {:s} with multiple labels."
                        .format(segm[0]))
            else:
                donoridd[segm[0]] = fld[0]
    return donoridd

def convert_id_list_str(id_list_str, iddict, fullidd, oufh = sys.stderr):
    newids = []
    fld = id_list_str.split(",")
    for id in fld:
        try:
            new_id = iddict[id]
        except KeyError:
            new_id = "NO_DONOR_{:s}".format(id)
        else:
            try:
                new_id = fullidd[new_id]
            except KeyError:
                new_id = "NO_GENOTYPE_{:s}".format(id)
        newids.append(new_id)

    oustr = "{:s}".format(newids[0])
    oufh.write(oustr + '\n')
    if len(newids) > 1:
        for new_id in newids[1:]:
            oustr += ',' + new_id
            oufh.write(new_id + '\n')
    return oustr

if __name__ == '__main__':
    nargs = len(sys.argv)

    if nargs < 3 or nargs > 5:
        sys.exit("usage: {} <conversion table> <comma separated list of ids> <file with full donor ids> [<output file>]".format(sys.argv[0]))

    fnam_conversion_table = sys.argv[1]
    id_list_str = sys.argv[2]
    fnam_full_donor_ids = sys.argv[3]
    if nargs > 4:
        oufh = open(sys.argv[4], 'w')
    else:
        oufh = sys.stderr

    iddict = load_conversion_table(fnam_conversion_table)
    donoridd = load_full_vcf_ids(fnam_full_donor_ids)
    oustr = convert_id_list_str(id_list_str, iddict, donoridd, oufh = oufh)

    sys.stdout.write(oustr + '\n')

    sys.exit(0)
