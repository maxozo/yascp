#!/usr/bin/env python


__date__ = '2022-10-20'
__author__ = 'M.Ozols'
__version__ = '0.0.1'

import argparse
import pandas as pd
import re
ibd_val=0.3
def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Joins all the outputs from the GT match in an expected input format for subsequent relatedness estimation.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
       '-ph', '--pihat_match',
        action='store',
        dest='pihat_match',
        required=True,
        help='PiHAT values for match'
    )

    parser.add_argument(
       '-c', '--cond',
        action='store',
        dest='condition',
        required=True,
        help='colname for condition'
    )
    
    parser.add_argument(
       '-e', '--expected',
        action='store',
        dest='expected',
        required=True,
        help='expected sample IDs'
    )
    
    
    parser.add_argument(
       '-md', '--metadata',
        action='store',
        dest='metadata',
        required=False,
        default=None,
        help='Sample metadata that mey contain relatedness info in DRAW DATE'
    )
        
    parser.add_argument(
       '-m', '--mapping_file',
        action='store',
        dest='mapping_file',
        required=False,
        default=None,
        help='Mapping between genotype/phenotype'
    )
    
    parser.add_argument(
       '-id', '--pool_id',
        action='store',
        dest='pool_id',
        required=False,
        default=None,
        help='Mapping between genotype/phenotype'
    )
    
    parser.add_argument(
       '-wpi', '--within_pool_ibd',
        action='store',
        dest='within_pool_ibd',
        required=False,
        default=None,
        help='Infered genotype relationships'
    )   

    parser.add_argument(
       '-mt', '--match_table',
        action='store',
        dest='match_table',
        required=True,
        help='match_table'
    )

    options = parser.parse_args()
    Genome_PiHAT = pd.read_csv(options.pihat_match,sep='\s+')
    GT_Match_Table = pd.read_csv(options.match_table,sep=',')
    # Name = options.match_table.replace('stats_','').replace('_gt_donor_assignments.csv','')
    Name = options.pool_id
    Condition_Column = options.condition
    expected_ids = options.expected.split(',')
    mapping_file = options.mapping_file
    wpi = options.within_pool_ibd
    if Condition_Column=='Expected' and wpi:    
        
        Genome_PiHAT_wpi = pd.read_csv(wpi,sep='\s+')
        Related_Pairs = Genome_PiHAT_wpi[Genome_PiHAT_wpi['PI_HAT']>ibd_val]
        all_unique_related_individuals = {}
        for id,row1 in Related_Pairs.iterrows():
            rid1 = row1['IID1']
            rid2 = row1['IID2']

            try:
                all_unique_related_individuals[rid1].append(rid2)
            except:
                all_unique_related_individuals[rid1]=[]
                all_unique_related_individuals[rid1].append(rid2)
                all_unique_related_individuals[rid1].append(rid1)
                
            try:
                all_unique_related_individuals[rid2].append(rid2)
            except:
                all_unique_related_individuals[rid2]=[]
                all_unique_related_individuals[rid2].append(rid2)
                all_unique_related_individuals[rid2].append(rid1)
    else:
        all_unique_related_individuals= {}
    
    if mapping_file:
        bridge = pd.read_csv(mapping_file,sep='\t')
        all_maped_samples = []
        for s1 in expected_ids:
            # there are no dublicates atm for the S2 ids, but there are replicates for the genotype ids, this however not conceir us.
            
            mapping = bridge[bridge['s00046_id']==s1]
            if len(mapping)==0:
                s2 = re.sub('^0*','',s1)
                mapping = bridge[bridge['s00046_id']==s2]
            if len(mapping)==0:
                s2 = s1.split('_')[0]
                mapping = bridge[bridge['s00046_id']==s2]
             
            if (len(mapping['oragene_id'])>0):
                for replacement in mapping['oragene_id']:
                    all_maped_samples.append({'original':s1,'replacement':replacement})
            else:
                replacement = s1
                all_maped_samples.append({'original':s1,'replacement':replacement})
        all_maped_samples2 = pd.DataFrame(all_maped_samples)
        all_maped_samples2 = all_maped_samples2.drop_duplicates()
        # all_maped_samples2=all_maped_samples2.set_index('original')
    else:
        all_maped_samples2 = expected_ids
        

    GT_Match_Table[f"PiHat: {Condition_Column}"]=None
    if Condition_Column=='Expected' and wpi:    
        GT_Match_Table[f"Infered Relatednes (PiHAT>{ibd_val})"]=None
    for i,row1 in GT_Match_Table.iterrows():
        donor_gt = row1['donor_gt']
        donor_querry = row1['donor_query']

        set1 =Genome_PiHAT[Genome_PiHAT['IID1']==donor_querry]
        set2 =Genome_PiHAT[Genome_PiHAT['IID2']==donor_querry]
        combo=pd.concat([set1,set2])

        set2_1 =combo[combo['IID1']==donor_gt]
        set2_2 =combo[combo['IID2']==donor_gt]
        combo2=pd.concat([set2_1,set2_2])
        
        if Condition_Column=='Expected' and wpi:    
            try:
               family_members =  all_unique_related_individuals[donor_querry]
               GT_Match_Table.loc[i,f"Infered Relatednes (PiHAT>{ibd_val})"]=';'.join(family_members)
            except:
                _=''
                
        try:
            GT_Match_Table.loc[i,f"PiHat: {Condition_Column}"]=combo2['PI_HAT'].values[0]
        except:
            _=''


                

    GT_Match_Table.to_csv(f"PiHAT_Stats_File_{Name}.csv",index=False,sep=',')


    if Condition_Column=='Expected':
        # Here we produce a new file that indicates all the best matches for the expected samples and their associated PiHat values.
        
        if (options.metadata):
            # Name = 'CRD_CMB12922404'
            Relationship_Clusters = {}
            Extra_Metadata_File = pd.read_csv(options.metadata,sep='\t')
            try:
                Pool_Metadata = Extra_Metadata_File[Extra_Metadata_File.experiment_id.str.contains(Name)]
                Pool_Metadata = Pool_Metadata.dropna(subset=['FAMILY','DRAW_DATE'])
                list_families = Pool_Metadata['FAMILY'].fillna('')+' '+Pool_Metadata['DRAW_DATE']
                
                list_families.index= Pool_Metadata['donor']
                list_families = list_families.dropna()
                Unique_Families = set(list_families)
            
                for family_id in Unique_Families:
                    Related_individuals = set(list_families[list_families == family_id].index.values)
                    for renated_indiv in Related_individuals:
                        Relationship_Clusters[renated_indiv]=Related_individuals
                        

                if len(list_families.dropna())>0:
                    perform_family_check = True
                else:
                    perform_family_check = False
            except:
                perform_family_check = False
                
        all_expected_id_best_match_PiHat = []
        for i,s1 in all_maped_samples2.iterrows():
            if (options.metadata):
                
                if perform_family_check:
                    try:
                        Reported_Relatives = Relationship_Clusters[s1.original]
                        
                        # Reported_Relatives.remove(s1.original)
                        relatedness = ';'.join(Reported_Relatives)
                    except:
                        relatedness = 'No Reports'
                else:
                    relatedness = 'No Reports'
            else:
                relatedness = 'No Reports'
 
            
            set1 =Genome_PiHAT[Genome_PiHAT['IID1'].str.contains(s1.replacement)]
            set2 =Genome_PiHAT[Genome_PiHAT['IID2'].str.contains(s1.replacement)]
            combo=pd.concat([set1,set2])
            try:
                max_pihat = combo['PI_HAT'].values.max()
                matching_donor_row = combo[combo['PI_HAT'] == max_pihat]
                if (matching_donor_row['IID1'].values[0]==s1.replacement):
                    donor_IID = matching_donor_row['IID2'].values[0]
                else:
                    donor_IID = matching_donor_row['IID1'].values[0]
                
            except:
                if (s1.original == s1.replacement):
                    donor_IID = 'mapping file seems to not contain this sample'
                else:
                    donor_IID = 'likely that sample mapping indicated is not in VCF file'
                max_pihat = None
            try:
                donor_nr =int(donor_IID.replace('donor',''))
            except:
                donor_nr = -1
            if Condition_Column=='Expected' and wpi:    
                try:
                    family_members =  all_unique_related_individuals[donor_IID]
                    Family_relatedness= ';'.join(family_members)
                except:
                    Family_relatedness= None
            else:
               Family_relatedness= None 
            all_expected_id_best_match_PiHat.append({'expected_donor_id':s1.original,'Best_Matched_donor_IID':donor_IID,'PiHAT value':max_pihat, 'mapping_id':s1.replacement,'donor_nr':donor_nr, 'reported relatedness':relatedness,f'infered genotype relatedness (PiHAT>{ibd_val})':Family_relatedness})
                
        all_expected_data_PiHAT = pd.DataFrame(all_expected_id_best_match_PiHat)
        all_expected_data_PiHAT= all_expected_data_PiHAT.sort_values(by=['donor_nr','PiHAT value'])
        all_expected_data_PiHAT=all_expected_data_PiHAT.drop(columns=['donor_nr'])
        all_expected_data_PiHAT.to_csv(f'Max_PiHAT_For_Expected_{Name}.tsv',sep='\t',index=False)
        


if __name__ == '__main__':
    main()
