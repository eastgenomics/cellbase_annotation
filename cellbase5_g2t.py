from pandas.core.frame import DataFrame
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import pandas as pd
from pandas.core.frame import DataFrame

from host import host_address

## Query cellbase5 database 
custom_config = {'rest': {'hosts': [host_address]}, 'version': 'v5', 'species': 'hsapiens'}
cc = ConfigClient(custom_config)
cbc = CellBaseClient(cc)
cbc.show_configuration()['version']
gc = cbc.get_gene_client() #select gene clients

# Read in the HGNC list that currently exists 
HGNC_df = pd.read_csv('HGNC_210902.tsv', sep='\t')
# this list will contain all the HGNC_ids and relevant MANE RefSeq transcripts
all_genes = []


for row in range(len(HGNC_df)):
    HGNC_id = HGNC_df.iloc[row, 0] #HGNC ID is in 0th column
    ensemble_id = HGNC_df.iloc[row, 11] #Ensembl gene ID is in 12th column
    # some genes do not have ensemble id so skip these
    if pd.isna(ensemble_id):
        print(HGNC_id + " does not have ensemble id to gene")
        pass
    else:
        # make dict to keep info on what refseq is for ensemble gene in cellbase
        gene_dict = {}
        gene_dict['HGNC_ID'] = HGNC_id
        # get all info and select transcript info
        data = gc.get_info(ensemble_id)
        # some ensemble gene id not present in cellbase so skip these
        if not data['responses'][0]['results']:
            print(ensemble_id + " does not exist in cellbase5") 
        else:
            # how many transcripts are there?
            txs_num = len(data['responses'][0]['results'][000]['transcripts'])
            # for loop the ensemble transcripts to find what the refseq mane 
            # transcripts are for the ensemble gene
            for transcript in range(txs_num):
                dicts = list(data['responses'][0]['results'][0]['transcripts'][transcript]['xrefs'])
                mane_dict = [item for item in dicts if item["dbName"] == "mane_select_refseq"]
                if not mane_dict:
                    # some ensemble transcripts do not have refseq mane transcript
                    # for these, we skip to the next transcript
                    continue
                else:
                    # refseq mane id is selected and added to the dict
                    gene_dict['MANE_RefSeqID'] = mane_dict[0]['id']
            # append the gene_dict to the list of all HGNC ids
            all_genes.append(gene_dict) 



# Save all hgnc ids and relevent refseq mane transcripts to dataframe 
df = pd.DataFrame.from_dict(all_genes)
df_noNaN = df[pd.notnull(df['MANE_RefSeqID'])]
# rest index to allow merging of columns later
df_noNaN = df_noNaN.reset_index(drop=True) 

#make two lists for canonical and clinical
clinical_list = ['clinical_transcript'] * len(df_noNaN)
canonical_list = ['canonical'] * len(df_noNaN)

#change into dataframe columns
clinical_list = pd.DataFrame({'clinical_list': clinical_list})
canonical_list = pd.DataFrame({'canonical_list': canonical_list})

#merge dataframe and columnc
dataframe = pd.concat([df_noNaN, clinical_list,canonical_list], axis=1)

# save output
dataframe.to_csv("g2t_210907_b38.tsv", sep="\t", header=False, index=False)

df.to_csv("g2t_210907_b38_all.tsv", sep="\t", header=False, index=False)