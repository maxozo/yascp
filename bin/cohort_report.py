#!/usr/bin/env python3
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(
    description="""
        Combines all the Raw cellranger metadata to be passed in the h5ad
        """
)

parser.add_argument(
    '-d',
    action='store',
    dest='path',
    required=False,
    default='./results',
    help='Path')

options = parser.parse_args()
path = options.path

try:
    prefix=f'{path}/handover/' #this is used for the updating reports posthoc, -> for this disable the next line
    project_name = os.listdir(f"{prefix}/Summary_plots")[0]
    GT_MATCH = pd.read_csv(f"{path}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
    Donor_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Donor_Report.tsv",sep='\t')
    Tranch_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Tranche_Report.tsv",sep='\t')
    Extra_Metadata_Donors = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
except:
    prefix='.'
    project_name = os.listdir(f"{prefix}/Summary_plots")[0]
    GT_MATCH = pd.read_csv(f"{path}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
    Donor_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Donor_Report.tsv",sep='\t')
    Tranch_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Tranche_Report.tsv",sep='\t')
    Extra_Metadata_Donors = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
# Split Donor Report by cohort


Donor_Report2 = Donor_Report.set_index('Pool_ID.Donor_Id')
try:
    Donor_Report2.insert(0,'Vacutainer ID','NONE')
except:
    print('Vacutainer ID exists')
# here if we enable we can replace the expected number of samples if a mistake was made.
# nr ukbb samples is already in there.

# qc/Cardinal_45327_Jul_18_2022/work/62/aacd091ea7711fae264e884e122230/Summary_plots/Cardinal_45327_Jul_18_2022/Summary
# t2 = GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline')]
# GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline'),'Match Expected']=True
GT_MATCH_CONFIDENT = GT_MATCH[GT_MATCH['Match Expected'] == True]

def Generate_Report(GT_MATCH_CONFIDENT,pan):

    Total_Report = pd.DataFrame()               
    for i,row1 in GT_MATCH_CONFIDENT.iterrows():
        print(i)
        # print(row1)
        print(f"{row1['pool']}_{row1['donor_query']}")
        
        Matched_Donor_report = ''
        Matched_Donor_report = Donor_Report[Donor_Report['Pool_ID.Donor_Id']==f"{row1['pool']}_{row1['donor_query']}"]
        
        Tranch_stats = Tranch_Report[Tranch_Report['Pool id'] ==row1['pool']]
        for i,row12 in Matched_Donor_report.iterrows():
            print(i)
            print(row12)
        if ('celline' in row1['donor_gt']):
            print('celline')
            gt_match = row1['donor_gt'].split('___')[1].replace('celline_','')
        else:
            gt_match = row1['donor_gt']
        
        Extra_Metadata_Donor1 = Extra_Metadata_Donors[Extra_Metadata_Donors['experiment_id'].str.contains(f"{row1['pool']}__{gt_match}")]
        try:
            Matched_Donor_report.insert(2, "Vacutainer ID", gt_match)
        except:
            print('exists')
        Matched_Donor_report['Experiment ID'] = project_name
        Matched_Donor_report['Chromium channel number'] = Tranch_stats['Chromium channel number'].values[0]
        Matched_Donor_report['Date of sample sequencing'] = Tranch_stats['Date of sample sequencing'].values[0]
        cohort = row1['final_panel']
        try:
            live_cell_count = Extra_Metadata_Donor1['live_cell_count'].values[0]
            viability =Extra_Metadata_Donor1['viability'].values[0]
            SITE = Extra_Metadata_Donor1['SITE'].values[0]
            Amount = Extra_Metadata_Donor1['Amount'].values[0]
            RECIEVED = Extra_Metadata_Donor1['RECIEVED'].values[0]
        except:
            live_cell_count = 'NA'
            viability ='NA'
            SITE = 'NA'
            Amount = 'NA'
            RECIEVED = 'NA'
            
        try:
            Matched_Donor_report.insert(8, 'lab_live_cell_count',live_cell_count)
        except:
            print('exists')
            Matched_Donor_report['lab_live_cell_count'] =  live_cell_count

        try:
            Matched_Donor_report.insert(8, 'viability', viability)
        except:
            print('exists')
            Matched_Donor_report['viability'] =  viability

        try:
            Matched_Donor_report.insert(8, 'cohort', cohort)
        except:
            print('exists')
            Matched_Donor_report['cohort'] =  cohort

        try:
            Matched_Donor_report.insert(8, 'Match Expected', row1['Match Expected'])
        except:
            print('exists')
            Matched_Donor_report['Match Expected'] =  row1['Match Expected']

        try:
            Matched_Donor_report.insert(8, 'site', SITE)
        except:
            print('exists')
            Matched_Donor_report['site'] =  SITE

        try:
            Matched_Donor_report.insert(8, 'amount recieved', Amount)
        except:
            print('exists')
            Matched_Donor_report['amount recieved'] =  Amount

        Matched_Donor_report['Date sample received'] =  RECIEVED
        Matched_Donor_report['Donor id']=row1['donor_gt original'].split('_')[0]
        Total_Report = pd.concat([Total_Report,Matched_Donor_report])

    return Total_Report


for confident_panel in set(GT_MATCH['final_panel']):
    if (confident_panel!='NONE' and confident_panel!='GT_cell_lines'):
        print(confident_panel)
        pan=confident_panel.replace('GT_','')
        print(pan)
        GR_PANEL = GT_MATCH[GT_MATCH['final_panel'] == confident_panel]
        GT_CELLINE = GT_MATCH[GT_MATCH['final_panel'] == 'GT_cell_lines']
        for i,row1 in GT_CELLINE.iterrows():
            celline = row1['donor_gt original'].split('_')[1]
            if celline in row1['Good_ids expected']:
                GT_CELLINE.loc[i,'Match Expected']=True

        GT_PANEL_CELLINE = pd.concat([GR_PANEL,GT_CELLINE])

        Total_Report = Generate_Report(GT_PANEL_CELLINE,pan)
        Total_Report2 = Total_Report.set_index('Pool_ID.Donor_Id')
        Donor_Report2.loc[Total_Report2.index,'Vacutainer ID']=Total_Report2.loc[Total_Report2.index,'Vacutainer ID']
        Donor_Report2.loc[Total_Report2.index,'Chromium channel number']=Total_Report2.loc[Total_Report2.index,'Chromium channel number']

        try: 
            Donor_Report2.insert(12, "site",'NONE') 
            Donor_Report2.insert(12, "lab_live_cell_count",'NONE') 
            Donor_Report2.insert(12, "viability",'NONE') 
        except: 
            print('exists')

        Donor_Report2.loc[Total_Report2.index,'site']=Total_Report2.loc[Total_Report2.index,'site']
        Donor_Report2.loc[Total_Report2.index,'lab_live_cell_count']=Total_Report2.loc[Total_Report2.index,'lab_live_cell_count']
        Donor_Report2.loc[Total_Report2.index,'viability']=Total_Report2.loc[Total_Report2.index,'viability']
        Donor_Report2 = Donor_Report2.reset_index().set_index('Pool ID')

        GT_MATCH2 = GT_MATCH.set_index('pool')
        GT_MATCH2=GT_MATCH2[['Emergency_ids expected','Good_ids expected']]
        GT_MATCH2 = GT_MATCH2.drop_duplicates().fillna("")
        GT_MATCH2['All IDs expected'] = GT_MATCH2['Emergency_ids expected']+','+GT_MATCH2['Good_ids expected']
        Donor_Report2['All IDs expected']=' '
        for idx1 in GT_MATCH2.index:
            # print(idx1)
            Donor_Report2.loc[idx1,'All IDs expected']=GT_MATCH2.loc[idx1,'All IDs expected']
        Donor_Report2 = Donor_Report2.reset_index().set_index('Pool_ID.Donor_Id')

        if (pan=='UKBB'):
            # For UKB samples we only return the expected samples
            Total_Report = Total_Report[Total_Report['Match Expected']]

        if confident_panel == 'GT_ELGH':
            confident_panel_name = 'Cardinal ELGH'
        elif confident_panel == 'GT_UKBB':
            confident_panel_name = 'Cardinal UKB'
        Expected_Samples = Extra_Metadata_Donors[Extra_Metadata_Donors.cohort == confident_panel_name]
        Missing_Samples = set(Expected_Samples.donor)-set(Total_Report['Vacutainer ID'])
        try:
            os.mkdir(f'{prefix}/Summary_plots/{project_name}/Summary/{pan}_REPORT')
        except:
            print('Dir exists')
        if (len(Missing_Samples)>0):
            Extra_Metadata_Donors3 = Extra_Metadata_Donors.copy()
            Extra_Metadata_Donors3 = Extra_Metadata_Donors3.set_index('donor')
            Missing = Extra_Metadata_Donors3.loc[Missing_Samples]['experiment_id']
            Missing = Missing.drop_duplicates()
            Missing.to_csv(f'{prefix}/Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_Missing_{pan}_Donors.tsv',sep='\t')
        Total_Report.to_csv(f'{prefix}/Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_{pan}_Report.tsv',sep='\t',index=False)
     
Donor_Report2.to_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Donor_Report.tsv",sep='\t',index=True)
Donor_Report2.to_csv(f"{prefix}/Summary_plots/{project_name}/Summary/{project_name}_Donor_Report.tsv",sep='\t',index=True)

print('Done')