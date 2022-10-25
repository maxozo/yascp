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
# input_id = 'ELGH_VAL11650907'
# stats_CRD_CMB12979963_gt_donor_assignments.csv
GT_Assignments = pd.read_csv(f'stats_{input_id}_gt_donor_assignments.csv',sep=',')
GT_Assignments=GT_Assignments.set_index('donor_query')
import os
tranche = os.getcwd().split('/')[-4]
pool = input_id
GT_Assignments__sample_summary_txt = pd.read_csv(f'{input_id}.sample_summary.txt',sep='\t',header=None)
GT_Assignments__exp_sample_summary_txt = pd.read_csv(f'{input_id}__exp.sample_summary.txt',sep='\t',header=None)
donor_ids=pd.read_csv(f'donor_ids.tsv',sep='\t')
vcf = pd.read_csv(f'GT_donors.vireo.vcf',sep='XXXXXXXdummy')
input_table_file = pd.read_csv(args.input_file,sep='\t')
input_table_file =input_table_file.set_index('experiment_id')
slipt1 = vcf[:1].values[0,0].split('\t')
# detect the outlier z values and make this as unassigned.

# here we use a different approach to assign confidence score
# from scipy.stats import iqr
# iqr = iqr(GT_Assignments['z0'])
# q3 = np.percentile(GT_Assignments['z0'], 75) 
# q1 = np.percentile(GT_Assignments['z0'], 25) 
# outlier_thresh = q3+1.5*iqr
# outlier_thresh_neg = q1-1.5*iqr

# Unassigned = list(GT_Assignments[GT_Assignments['z0']<outlier_thresh_neg].index)
# Unassigned2 = list(GT_Assignments[GT_Assignments['z0']<50].index)
# Unassigned1 = list(GT_Assignments[GT_Assignments['z0']-GT_Assignments['z1']<=2].index)
# Unassigned= list(set(Unassigned2+Unassigned1))


Unassigned= []
All_expected_ids = input_table_file.loc[input_id,'donor_vcf_ids'].replace('\'','').split(',')
Good_ids=[]
Emergency_ids=[]

for id1 in All_expected_ids:
    if '-999-' in id1:
        Emergency_ids.append(id1)
    else:
        Good_ids.append(id1)
Good_ids = ",".join(Good_ids)
Emergency_ids = ",".join(Emergency_ids)
GT_Assignments['Emergency_ids expected']=Emergency_ids
GT_Assignments['Good_ids expected']=Good_ids
GT_Assignments['Match Expected']='False'
GT_Assignments['pool']=pool
GT_Assignments['tranche']=tranche
GT_Assignments['donor_gt original']=GT_Assignments['donor_gt']
D2 = pd.DataFrame(All_expected_ids,columns=['col1'])
for ix in GT_Assignments.index:
    # print(ix)
    
    #Here should add a a filter to estimate whether it is a good match.
    replacement = GT_Assignments.loc[ix,'donor_gt']

    expected = 'NA'
    poor_replacement =''
    if ix in Unassigned:
       poor_replacement = 'Poor_GT_score___'

    # The folowing will attempt to replcat the donor names as per phenotype file, as provided by mapping.
    # If it fails it will keep the genotype id.
    
    if (len(gp_ma)>0):
        try:
            # replacement = 'S2-046-01741'
            replacement = str(gp_ma.loc[replacement][0])
            print(replacement)
            # replacement = replacement.values[0]
            
            if len(D2[D2.col1.str.contains(replacement)])>0:
                GT_Assignments.loc[ix,'Match Expected']='True'
        except:
            try:
                # replacement='15001608190388_204238910153_R09C02'
                
                replacements = gp_ma.loc[replacement.split('_')[0]]
                replacement = 'No_mapping___'+replacement
                # replacement = replacement.values[0]
                for rep1 in replacements.iloc[:,0]:
                    
                    if len(D2[D2.col1.str.contains(rep1)])>0:
                        GT_Assignments.loc[ix,'Match Expected']='True'
                        replacement = rep1
            except:
                replacement = 'No_mapping___'+replacement
    else:
        _ = ''
    if 'THP1' in replacement:
        replacement = 'celline_THP1'
    if 'U937' in replacement:
        replacement = 'celline_U937'
    if poor_replacement == 'Poor_GT_score___':
        replacement = poor_replacement+replacement

    GT_Assignments__exp_sample_summary_txt.loc[:,0] = GT_Assignments__exp_sample_summary_txt.loc[:,0].replace(f'{input_id}__{ix}', f'{input_id}__{replacement}')
    GT_Assignments__sample_summary_txt.loc[:,1]=GT_Assignments__sample_summary_txt.loc[:,1].replace(ix, replacement)
    donor_ids['donor_id']=donor_ids['donor_id'].replace(ix, replacement)
    donor_ids['best_singlet']=donor_ids['best_singlet'].replace(ix, replacement)
    sp1 = donor_ids['best_doublet'].str.split(',').str[0].replace(ix, replacement)
    sp2 = donor_ids['best_doublet'].str.split(',').str[1].replace(ix, replacement)
    donor_ids['best_doublet'] =sp1+','+sp2
    slipt1=list(map(lambda x: x.replace(ix, replacement), slipt1))
    GT_Assignments.loc[ix,'donor_gt']=replacement


slipt1='\t'.join(slipt1)
vcf[:1].values[0,0]=slipt1
# Now that we have loaded them replace all the files with the correct donor ids.

GT_Assignments__sample_summary_txt.to_csv(f'GT_replace_{input_id}.sample_summary.txt',sep='\t',index=False)
GT_Assignments__exp_sample_summary_txt.to_csv(f'GT_replace_{input_id}__exp.sample_summary.txt',sep='\t',index=False)
GT_Assignments.to_csv(f'GT_replace_{input_id}_assignments.tsv',sep='\t')
donor_ids.to_csv(f'GT_replace_donor_ids.tsv',sep='\t',index=False)
vcf.to_csv(f'GT_replace_GT_donors.vireo.vcf',sep=',',index=False)

print("Done")