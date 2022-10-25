#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np

from scipy import stats

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-id", "--donor_id", metavar="donor_id", dest = "donor_id",
                    help="", type=str,
                    default=1)
parser.add_argument("-gp", "--genotype_phenotype_mapping", metavar="genotype_phenotype_mapping", dest = "genotype_phenotype_mapping",
                    help="", type=str,
                    default=None)

parser.add_argument("-if", "--input_file", metavar="input_file", dest = "input_file",
                    help="", type=str,
                    default=None)
                    
args = parser.parse_args()
if (args.genotype_phenotype_mapping):
    gp_ma = pd.read_csv(args.genotype_phenotype_mapping,sep='\t',index_col=0)
    gp_ma.index = gp_ma.index.astype(str)
else:
    gp_ma = pd.DataFrame()
input_id =args.donor_id
donors_in_vcf = pd.read_csv('donors_in_vcf.tsv',header=None)
# input_id = 'ELGH_VAL11650907'
# stats_CRD_CMB12979963_gt_donor_assignments.csv
# GT_Assignments = pd.read_csv(f'stats_{input_id}_gt_donor_assignments.csv',sep=',')
# GT_Assignments=GT_Assignments.set_index('donor_query')
import os
tranche = os.getcwd().split('/')[-4]
pool = input_id
GT_Assignments__sample_summary_txt = pd.read_csv(f'{input_id}.sample_summary.txt',sep='\t',header=None)
GT_Assignments__exp_sample_summary_txt = pd.read_csv(f'{input_id}__exp.sample_summary.txt',sep='\t',header=None)
donor_ids=pd.read_csv(f'donor_ids.tsv',sep='\t')
vcf = pd.read_csv(f'GT_donors.vireo.vcf.gz',sep='XXXXXXXdummy')
input_table_file = pd.read_csv(args.input_file,sep='\t')
input_table_file =input_table_file.set_index('experiment_id')
slipt1 = vcf[:1].values[0,0].split('\t')
GT_Assignments = pd.DataFrame()
# s2 = pd.Series(['donor3','donor5'])
GT_Assignments['ix'] = GT_Assignments__sample_summary_txt.loc[:,1]
GT_Assignments = GT_Assignments.reset_index(drop=True)
GT_Assignments['donor_gt'] = ''
donor_names_set = set(GT_Assignments['ix']).union(set(donors_in_vcf.iloc[:,0]))
All_Dons_pre = pd.DataFrame(donor_names_set,columns=['ix'])
GT_Assignments = All_Dons_pre
GT_Assignments['donor_gt'] = ''
All_Dons = All_Dons_pre[All_Dons_pre['ix'].str.contains('donor')]


int1 = 0

Taken_nrs = set(All_Dons['ix'].str.replace('donor','').astype(int))
for i,row1 in GT_Assignments.iterrows():
    don1 = row1['ix']
    if don1 in ['doublet','unassigned']:
        replacement = don1
    else:
        if ('donor' in don1):
            replacement =don1
        else:
            while int1 in Taken_nrs:
                int1+=1
            replacement = f"donor{int1}"
            int1+=1
    GT_Assignments.loc[i,'donor_gt'] = replacement
    GT_Assignments__exp_sample_summary_txt.loc[:,0] = GT_Assignments__exp_sample_summary_txt.loc[:,0].replace(f'{input_id}__{don1}', f'{input_id}__{replacement}')
    GT_Assignments__sample_summary_txt.loc[:,1]=GT_Assignments__sample_summary_txt.loc[:,1].replace(don1, replacement)
    donor_ids['donor_id']=donor_ids['donor_id'].replace(don1, replacement)
    donor_ids['best_singlet']=donor_ids['best_singlet'].replace(don1, replacement)
    sp1 = donor_ids['best_doublet'].str.split(',').str[0].replace(don1, replacement)
    sp2 = donor_ids['best_doublet'].str.split(',').str[1].replace(don1, replacement)
    donor_ids['best_doublet'] =sp1+','+sp2
    slipt1=list(map(lambda x: x.replace(don1, replacement), slipt1))
    


slipt1='\t'.join(slipt1)
vcf[:1].values[0,0]=slipt1
# Now that we have loaded them replace all the files with the correct donor ids.

GT_Assignments__sample_summary_txt.to_csv(f'GT_replace_{input_id}.sample_summary.txt',sep='\t',index=False)
GT_Assignments__exp_sample_summary_txt.to_csv(f'GT_replace_{input_id}__exp.sample_summary.txt',sep='\t',index=False)
GT_Assignments.to_csv(f'GT_replace_{input_id}_assignments.tsv',sep='\t')
donor_ids.to_csv(f'GT_replace_donor_ids.tsv',sep='\t',index=False)
# vcf.to_csv(f'GT_replace_GT_donors.vireo.vcf',sep=',',index=False)
GT_Assignments.to_csv(f'replacement_assignments.tsv',sep=' ',index=False,header=False)
print("Done")