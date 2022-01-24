import argparse
import datetime
import pandas as pd
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import sys
import time


def parse_args():
    """Reads on argument passed at the cmd line
    Returns:
        args: args held in an object
    """

    parser = argparse.ArgumentParser(
        description='HGNC tsv containing columns HGNC and ensembl gene ID')
    parser.add_argument(
        '-f', '--file',
        required=True,
        help='HGNC tsv'
        )

    args = parser.parse_args()

    return args


def query_cellbase(gc, HGNC_df, HGNC_missing_ensemblID,
                    ensemblID_not_in_cellbase,
                    ensemblID_has_no_maneselect_refseq, all_genes):
    """This function queries cellbase for a given gene from the HGNC
    table. As there is not get HGNC ID function for cellbase, the
    ensembl gene id will be queried instead. The output is the equivalent
    mane refseq transcript for the gene and three lists for HGNC id with
    no ensembl ID, ensembl gene ID not in cellbase and ensembl gene ID
    not assigned to a mane refseq transcript.
    Inputs:
        gc: cellbase gene client that is queried with a ensembl gene id
        HGNC_df: dataframe containing HGNC id and ensembl gene id column
        HGNC_missing_ensemblID: empty list for HGNC ID
                                that dont have an ensembl gene
        ensemblID_not_in_cellbase: empty list for ensembl gene
                                    id not in cellbase
        ensemblID_has_no_maneselect_refseq: empty list for ensembl genes
                                        that dont have maneselect_refseq
        all_genes: empty list that will hold HGNC id and
                    their maneselect refseq transcript

    Returns:
        all_genes: list holding dictionaries of HGNC id and
                    their maneselect refseq transcript
        HGNC_missing_ensemblID: list for HGNC ID that dont
                                have an ensembl gene
        ensemblID_not_in_cellbase: list for ensembl gene id
                                    not in cellbase
        ensemblID_has_no_maneselect_refseq: list for ensembl genes that
                                            dont have maneselect_refseq
    """

    for row in range(len(HGNC_df)):
        # HGNC ID is in 1st column
        HGNC_id = HGNC_df.loc[row, 'HGNC ID']
        # Ensembl gene ID is in 12th column
        ensembl_id = HGNC_df.loc[row, 'Ensembl gene ID']

        # some genes do not have ensembl id so skip these
        if pd.isna(ensembl_id):
            print(HGNC_id + " does not have ensembl id to gene")
            HGNC_missing_ensemblID.append(HGNC_id)
            continue
        # make dict to keep info on what refseq is for ensembl gene in cellbase
        gene_dict = {}
        gene_dict['HGNC_ID'] = HGNC_id

        # get all info and select transcript info
        # since it will query from cellbase it may fail, its best
        # to raise an exception to catch it
        try:
            data = gc.get_info(ensembl_id)
        except:
            print(f"Unable to get {ensembl_id} info from cellbase cleint")

        # some ensembl gene id not present in cellbase so skip these
        if not data['responses'][0]['results']:
            print(f"{ensembl_id} does not exist in cellbase5")
            ensemblID_not_in_cellbase.append(ensembl_id)
        # how many transcripts are there?
        txs_num = len(data['responses'][0]['results'][000]['transcripts'])
        # for loop the ensembl transcripts to find what the refseq mane
        # transcripts are for the ensembl gene
        for transcript in range(txs_num):
            dicts = list(data['responses'][0]['results'][0]['transcripts'][transcript]['xrefs'])
            mane_dict = [item for item in dicts if item["dbName"] == "mane_select_refseq"]
            if not mane_dict:
                # some ensembl transcripts do not have refseq mane transcript
                # for these, we skip to the next transcript
                ensemblID_has_no_maneselect_refseq.append(ensembl_id)
                continue
            else:
                # refseq mane id is selected and added to the dict
                gene_dict['MANE_RefSeqID'] = mane_dict[0]['id']
        # append the gene_dict to the list of all HGNC ids
        all_genes.append(gene_dict)

    return (all_genes,
            HGNC_missing_ensemblID,
            ensemblID_not_in_cellbase,
            ensemblID_has_no_maneselect_refseq)


def main():
    """
    Main function to set up cellbase version 5, put the HGNC dataframe
    through the query_cellbase function and save all outuputs.
    """
    # time taken for script to run
    begin_time = datetime.datetime.now()

    # print all outputs to an output job file
    # named tuple contains local time
    named_tuple = time.localtime()
    today_datetime = time.strftime("%Y%m%d_%H%M%S", named_tuple)
    sys.stdout = open(today_datetime + "_output_log.txt", "w")

    # Query cellbase5 database
    custom_config = {'rest': {'hosts': [
                                'https://ws.zettagenomics.com/cellbase'
                                ]},
                                'version': 'v5', 'species': 'hsapiens'}
    customconfigclient = ConfigClient(custom_config)
    cbc = CellBaseClient(customconfigclient)
    cbc.show_configuration()['version']
    gc = cbc.get_gene_client()  # select gene clients

    # read in the input file hgnc tsv from the cmd line
    args = parse_args()
    print(
        "================      Input files     ================\n",
        f"HGNC table input file => {args.file}"
        )
    HGNC_df = pd.read_csv(args.file, sep='\t')
    all_genes = []

    # list of things that are missing for certain HGNC/ensembl IDs
    HGNC_missing_ensemblID = []
    ensemblID_not_in_cellbase = []
    ensemblID_has_no_maneselect_refseq = []

    # input the lists and get out the filled out lists
    all_genes, HGNC_missing_ensemblID, ensemblID_not_in_cellbase, \
    ensemblID_has_no_maneselect_refseq = query_cellbase(
        gc, HGNC_df, HGNC_missing_ensemblID, ensemblID_not_in_cellbase,
        ensemblID_has_no_maneselect_refseq, all_genes)

    # Save all hgnc ids and relevent refseq mane transcripts to dataframe
    df = pd.DataFrame.from_dict(all_genes)
    # some of df contains genes with no manerefseq ID,
    # they have NaNs so remove it
    df_noNaN = df[pd.notnull(df['MANE_RefSeqID'])]
    # rest index to allow merging of columns later
    df_noNaN = df_noNaN.reset_index(drop=True)

    # Add two lists for canonical and clinical
    df_noNaN['clinical_list'] = 'clinical_transcript'
    df_noNaN['canonical_list'] = 'canonical'

    # save outputs
    print("================      Output files    ================")

    df_noNaN_outputname = today_datetime + '_g2t_b38.tsv'
    df_outputname = today_datetime + '_g2t_b38_all.tsv'
    df_noNaN.to_csv(df_noNaN_outputname, sep="\t", header=False, index=False)
    df.to_csv(df_outputname, sep="\t", header=False, index=False)

    print(
        "Table with HGNC IDs with mane_select refseq transcript ",
        f"=> {df_noNaN_outputname}"
            )
    print(
        "Table with all HGNC IDs with and without mane_select ",
        f"refseq transcript => {df_outputname}"
            )

    # save list that wont be in the final file
    # change into dataframe columns
    HGNC_missing_ensemblID = pd.DataFrame({
        'HGNC_missing_ensemblID': HGNC_missing_ensemblID
        })
    ensemblID_not_in_cellbase = pd.DataFrame({
        'ensemblID_not_in_cellbase': ensemblID_not_in_cellbase
        })
    ensemblID_has_no_maneselect_refseq = pd.DataFrame({
        'ensembl_gene_has_no_maneselect_refseq': ensemblID_has_no_maneselect_refseq
        })
    # merge dataframe and column
    missinginfo_df = pd.concat(
        [HGNC_missing_ensemblID, ensemblID_not_in_cellbase,
        ensemblID_has_no_maneselect_refseq],
        axis=1
        )

    missinginfo_df_outputname = today_datetime + '_g2t_b38_missing_info.tsv'
    missinginfo_df.to_csv(missinginfo_df_outputname,
                        sep="\t", header=True, index=False)

    print(
        "Table with HGNC IDs missing ensembl gene IDs ",
        "ensembl gene IDs not in cellbase and ",
        f"=> {missinginfo_df_outputname}\n"
        )

    print(f"Execution time {(datetime.datetime.now() - begin_time)}")


if __name__ == "__main__":

    main()
