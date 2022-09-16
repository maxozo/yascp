#!/usr/bin/env python3

## read donor_ids from a file and compare them to a list of donor ids obtained from VCF header
## output a table with 2 columns <donor_id>,<present_in_vcf [Y/N]>


import sys

def read_list_file(fnam):
    items = []
    with open(fnam, 'r') as infh:
        for lin in infh:
            fld = lin.split()
            id = fld[0]
            if len(id) > 0:
                items.append(id)
    return items

if __name__ == '__main__':
    nargs = len(sys.argv)

    if nargs != 4:
        sys.exit("usage: {:s} <list_of_donor_ids_in_vcf> <list_of_donor_ids_in_pool> <output file [TSV]>"
            .format(sys.argv[0]))

    vcf_donor_file = sys.argv[1]
    pooled_donor_file = sys.argv[2]
    fnam_output = sys.argv[3]

    vcf_donors = read_list_file(vcf_donor_file)
    pool_donors = read_list_file(pooled_donor_file)

    oufh = open(fnam_output, 'w')
    oufh.write("donor_id\tis_in_VCF\n")
    for donor_id in pool_donors:
        if donor_id in vcf_donors:
            is_present = 'Y'
        else:
            is_present = 'N'
        oufh.write("{0:s}\t{1:s}\n".format(donor_id, is_present))
    oufh.close()

    sys.exit(0)
