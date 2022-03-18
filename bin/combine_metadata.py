#!/usr/bin/env python3

__date__ = '2021-04-11'
__version__ = '0.0.1'

import argparse
import pandas as pd



def main(files):
    dataframe =[]
    for file in files.split(","):
        experiment_id = file.split('---')[0]

        Data= pd.read_csv(file)
        data={}
        data['experiment_id']=experiment_id
        for col1 in Data.columns:
            col1_2 = col1.replace(' ','_')
            try:
                col1_val = Data[col1].values[0].replace(',','')
            except:
                col1_val = Data[col1].values[0]
            data[col1_2]=col1_val
        dataframe.append(data)
        Data2 = pd.DataFrame(dataframe)
        Data2.to_csv("full_metadata.tsv",sep='\t',index=False)



if __name__ == '__main__':


    parser = argparse.ArgumentParser(
        description="""
            Combines all the Raw cellranger metadata to be passed in the h5ad
            """
    )

    parser.add_argument(
        '-d',
        action='store',
        dest='paths',
        required=True,
        help='Metadata Paths')

    options = parser.parse_args()
    files = options.paths
    main(files)
