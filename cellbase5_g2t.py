from pandas.core.frame import DataFrame
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import pandas as pd
from pandas.core.frame import DataFrame
import argparse
from datetime import date

from host import host_address


def parse_args():
    """Reads on argument passed at the cmd line
    Returns:
        args: args held in an object
    """
    
    parser = argparse.ArgumentParser(
        description = 'HGNC tsv containg columns HGNC and ensemble gene ID')
    parser.add_argument(
        '-f', '--file',
        required=True,
        help='HGNC tsv'
        )

    args = parser.parse_args()

    return args

def query_cellbase(gc, HGNC_df, HGNC_missing_ensembleID, ensembleID_not_in_cellbase, ensembleID_has_no_maneselect_refseq, all_genes):
    """
    Inputs:
        gc: cellbase gene client that is queried with a ensemble gene id
        HGNC_df: dataframe containing HGNC id and ensemble gene id column 
        HGNC_missing_ensembleID: empty list for HGNC ID that dont have an ensemble gene
        ensembleID_not_in_cellbase: empty list for ensemble gene id not in cellbase
        ensembleID_has_no_maneselect_refseq: empty list for ensemble genes that dont have maneselect_refseq
        all_genes: empty list that will hold HGNC id and their maneselect refseq transcript
    
    Returns:
        all_genes: list holding dictionaries of HGNC id and their maneselect refseq transcript
        HGNC_missing_ensembleID: list for HGNC ID that dont have an ensemble gene
        ensembleID_not_in_cellbase: list for ensemble gene id not in cellbase
        ensembleID_has_no_maneselect_refseq: list for ensemble genes that dont have maneselect_refseq
    """
    
    for row in range(len(HGNC_df)):
        HGNC_id = HGNC_df.loc[row, 'HGNC ID'] #HGNC ID is in 1st column
        ensemble_id = HGNC_df.loc[row, 'Ensembl gene ID'] #Ensembl gene ID is in 12th column
        # some genes do not have ensemble id so skip these
        if pd.isna(ensemble_id):
            print(HGNC_id + " does not have ensemble id to gene")
            HGNC_missing_ensembleID.append(HGNC_id)
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
                ensembleID_not_in_cellbase.append(ensemble_id)
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
                        ensembleID_has_no_maneselect_refseq.append(ensemble_id)
                        continue
                    else:
                        # refseq mane id is selected and added to the dict
                        gene_dict['MANE_RefSeqID'] = mane_dict[0]['id']
                # append the gene_dict to the list of all HGNC ids
                all_genes.append(gene_dict) 
    
    return all_genes, HGNC_missing_ensembleID, ensembleID_not_in_cellbase, ensembleID_has_no_maneselect_refseq



def main():
    """
    Main function to set up cellbase version 5, put the HGNC dataframe
    through the query_cellbase function and save all outuputs.
    """
    ## Query cellbase5 database 
    custom_config = {'rest': {'hosts': [host_address]}, 
                    'version': 'v5', 'species': 'hsapiens'}
    cc = ConfigClient(custom_config)
    cbc = CellBaseClient(cc)
    cbc.show_configuration()['version']
    gc = cbc.get_gene_client() #select gene clients

    # read in the hgnc tsv from the cmd line
    args = parse_args()
    input_filename = args.file

    HGNC_df = pd.read_csv(input_filename, sep='\t')
    all_genes = []

    # list of things that are missing for certain HGNC/ensemble IDs
    HGNC_missing_ensembleID = []
    ensembleID_not_in_cellbase = []
    ensembleID_has_no_maneselect_refseq = []
    
    # input the lists and get out the 
    all_genes, HGNC_missing_ensembleID, ensembleID_not_in_cellbase, ensembleID_has_no_maneselect_refseq = query_cellbase(
        gc, HGNC_df, HGNC_missing_ensembleID, ensembleID_not_in_cellbase, 
        ensembleID_has_no_maneselect_refseq, all_genes)

    # Save all hgnc ids and relevent refseq mane transcripts to dataframe 
    df = pd.DataFrame.from_dict(all_genes)
    # some of df contains genes with no manerefseq ID, they have NaNs so remove it
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
    complete_dataframe = pd.concat([df_noNaN, clinical_list,canonical_list], axis=1)

    ##### save outputs
    today = date.today()

    complete_dataframe_outputname = today.strftime("%Y%m%d") +'_g2t_b38.tsv'
    df_outputname = today.strftime("%Y%m%d") +'_g2t_b38_all.tsv'
    complete_dataframe.to_csv(complete_dataframe_outputname, sep="\t", header=False, index=False)
    df.to_csv(df_outputname, sep="\t", header=False, index=False)

    #save list that wont be in the final file
    #change into dataframe columns
    HGNC_missing_ensembleID = pd.DataFrame({'HGNC_missing_ensembleID': HGNC_missing_ensembleID})
    ensembleID_not_in_cellbase = pd.DataFrame({'ensembleID_not_in_cellbase': ensembleID_not_in_cellbase})
    ensembleID_has_no_maneselect_refseq = pd.DataFrame({'ensemble_gene_has_no_maneselect_refseq': ensembleID_has_no_maneselect_refseq})
    #merge dataframe and column
    missinginfo_dataframe = pd.concat([HGNC_missing_ensembleID, ensembleID_not_in_cellbase,ensembleID_has_no_maneselect_refseq], axis=1)

    missinginfo_dataframe_outputname = today.strftime("%Y%m%d") +'_g2t_b38_missing_info.tsv'
    missinginfo_dataframe.to_csv(missinginfo_dataframe_outputname, sep="\t", header=True, index=False)

if __name__ == "__main__":

    main()