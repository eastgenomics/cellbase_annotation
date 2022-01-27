import argparse
from datetime import date
import pandas as pd
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient
import urllib.request, json

from host import host_address

def extract_exons_info(exon_dict, gene_dict):
    gene_dict['chr'] = exon_dict["chromosome"]
    gene_dict['exon_start'] = exon_dict["genomicCodingStart"]
    gene_dict['exon_end'] = exon_dict["genomicCodingEnd"]
    gene_dict['exonNumber'] = exon_dict["exonNumber"]
    return gene_dict

customconfigclient = ConfigClient('config.json')
cbc = CellBaseClient(customconfigclient)
cbc.show_configuration()
gc = cbc.get_gene_client()  # select gene clients

HGNC_df = pd.read_csv('HGNC_210902_two.tsv', sep='\t')
all_genes = []
# list of things that are missing for certain HGNC/ensembl IDs
HGNC_missing_ensemblID = []
ensemblID_not_in_cellbase = []
ensembleID_with_no_transcript = []
transcript_with_no_exon = []
transcript_with_no_coding_exon = []
df = pd.DataFrame()

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
    # if error arises -> print another error messages
    try:
        data = gc.get_info(ensembl_id)
    except:
        print("Killed - Client lost connection with Cellbase")

    # some ensembl gene id not present in cellbase so skip these
    if not data['responses'][0]['results']:
        print(f"{ensembl_id} does not exist in cellbase5")
        ensemblID_not_in_cellbase.append(ensembl_id)
    # how many transcripts are there?
    txs_num = len(data['responses'][0]['results'][000]['transcripts'])
    print(txs_num)

    # Add the gene symbol
    gene_dict['gene_symbol'] = data['responses'][0]['results'][0]['name']

    # some genes dont have transcripts so skip these and note these
    if txs_num == 0:
        ensembleID_with_no_transcript.append(ensembl_id)
        print('ensembleID_with_no_transcript')
    else:
        # for loop the ensembl transcripts to extract each exon
        for transcript in range(txs_num):
            # Multiple transcripts means we make a list of dictionaries per gene
            list_txs_dict = []
            gene_dict['transcript_id'] = data['responses'][0]['results'][0]['transcripts'][transcript]['id']
            transcript_dict = list(data['responses'][0]['results'][0]['transcripts'][transcript])
            exons_num = len(data['responses'][0]['results'][0]['transcripts'][transcript]['exons'])
            # some transcripts dont have exons so skip these and note these
            txs = data['responses'][0]['results'][0]['transcripts'][transcript]['id']
            print(txs)
            if exons_num == 0:
                transcript_with_no_exon.append(txs)
                print('transcript_with_no_exon')

            # loop through each exon for a transcript and make a long list of dict
            for exon in range(exons_num):
                exon_dict = data['responses'][0]['results'][0]['transcripts'][transcript]['exons'][exon]
                # phase = -1 means that the whole exon is not a coding exon
                # if its not a coding exon, dont include and note this
                if exon_dict["phase"] != -1:
                    gene_dict = extract_exons_info(exon_dict, gene_dict)
                else:
                    ex = data['responses'][0]['results'][0]['transcripts'][transcript]['exons'][exon]['id']
                    transcript_with_no_coding_exon.append([txs, ex])
                    print('transcript_with_no_coding_exon')
                #print(gene_dict)
                list_txs_dict.append(gene_dict)
                print(list_txs_dict)

            transcript_table = pd.DataFrame(list_txs_dict)
            df = df.append(transcript_table)


print("AAAAAAAAAAAHHHHHHHHHHHH")
print(df)
df[["exon_start", "exon_end", "exonNumber"]] = df[["exon_start", "exon_end", "exonNumber"]].fillna(0.0).astype(int)
#df[["exon_start", "exon_end", "exonNumber"]] = df[["exon_start", "exon_end", "exonNumber"]].astype('int64')
df.to_csv("first_exons.tsv", sep="\t", header=True, index=False)

#merge missing info

missinginfo_df = pd.DataFrame([HGNC_missing_ensemblID, ensemblID_not_in_cellbase,
    ensembleID_with_no_transcript, transcript_with_no_exon,
    transcript_with_no_coding_exon]
    )
missinginfo_df = missinginfo_df.transpose() #To Transpose and make each rows as columns
missinginfo_df.columns=['HGNC_missing_ensemblID', 'ensemblID_not_in_cellbase',
    'ensembleID_with_no_transcript', 'transcript_with_no_exon',
    'transcript_with_no_coding_exon']

missinginfo_df_outputname = 'missing_info.tsv'
missinginfo_df.to_csv(missinginfo_df_outputname,
                    sep="\t", header=True, index=False)