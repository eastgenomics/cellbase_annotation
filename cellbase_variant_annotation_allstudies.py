## creating a popyulation frequency of myeloid snps to be used for somalier
## luckily cellbase has this incorporated into their database !! :)

# e.g https://ws.zettagenomics.com/cellbase/webservices/rest/v5/hsapiens/genomic/variant/search?assembly=grch38&count=false&region=1%3A1787324-1787444&limit=10&skip=0


from pandas.core.frame import DataFrame
import pandas as pd
import urllib.request, json
from host import host_address


# via url
chrX_regions = pd.read_csv("~/Documents/Projects/SexCheck/testing_data/myeloid/probes_GRCh38_myeo_v1.0_chrX.bed", 
                            sep = "\t", names=["chr", "start","end"])

variant_dataframe = pd.DataFrame()

for row in range(len(chrX_regions)):
    print(row)
    chr = chrX_regions.iloc[row, 0]
    # zetta doesnt like Chr when query chromosome so strip the chr
    chr = chr[3:]
    start = chrX_regions.iloc[row, 1]
    end = chrX_regions.iloc[row, 2]
    url_start = "/webservices/rest/v5/hsapiens/genomic/variant/search?assembly=grch38&count=false&region="
    url_middle = str(chr) + "%3A" + str(start) + "-" + str(end)
    url_end = "&limit=10&skip=0"
    url_address = host_address + url_start  + url_middle + url_end
    print(url_address)
    with urllib.request.urlopen(url_address) as url:
        data = json.loads(url.read().decode())
    
    # number of variants annotated
    variant_number = len(data['responses'][0]['results'])
    # loop through each of these and extract the populationFrequencines annotations
    # we can filter out study specific cases later, lets keep it all for now
    #print(variant_number)
    for variant in range(variant_number):
        #create a new variant dataframe and append to main
        populationFreq = data['responses'][0]['results'][variant]['annotation']['populationFrequencies']
        #there are different annotation studies and populations, extract them all
        annot_length = len(populationFreq)
        # prepare the lists and info to get out
        # get the chr, pos of this variant and * by the number of annotations
        var_chr = [data['responses'][0]['results'][variant]['chromosome']] * annot_length
        var_pos = [data['responses'][0]['results'][variant]['start']] * annot_length
        var_ID = ["."] * annot_length
        var_ref = []
        var_alt = []
        var_qual = ["."] * annot_length
        var_filter = ["PASS"] * annot_length
        var_study = []
        var_population = []
        var_AF = []
        # loop through each annotation info for this variant and extract all
        for annotation in range(annot_length):
            annot_study = populationFreq[annotation]
            var_study.append(annot_study["study"])
            var_population.append(annot_study["population"])
            var_ref.append(annot_study["refAllele"])
            var_alt.append(annot_study["altAllele"])
            var_AF.append(annot_study["altAlleleFreq"])
        #change into dataframe columns
        var_chr = pd.DataFrame({'#CHROM': var_chr})
        var_pos = pd.DataFrame({'POS': var_pos})
        var_ID = pd.DataFrame({'ID': var_ID})
        var_ref = pd.DataFrame({'REF': var_ref})
        var_alt = pd.DataFrame({'ALT': var_alt})
        var_qual = pd.DataFrame({'QUAL': var_qual})
        var_filter = pd.DataFrame({'FILTER': var_filter})
        var_study = pd.DataFrame({'study': var_study})
        var_population = pd.DataFrame({'population': var_population})
        var_AF = pd.DataFrame({'INFO': var_AF})
        #append all the annotations of that one variant
        dataframe = pd.concat([var_chr, var_pos, var_ID, var_ref,var_alt,var_qual,var_filter, var_study,var_population,var_AF], axis=1)
        #append annotations to the other/main variants table
        variant_dataframe = variant_dataframe.append(dataframe)
        


variant_dataframe['INFO'] = 'AF=' + variant_dataframe['INFO']
variant_dataframe['#CHROM'] = 'chr' + variant_dataframe['#CHROM'].astype(str)
variant_dataframe.to_csv("~/Documents/Projects/SexCheck/testing_data/myeloid/probes_myeloid_variant_popfreq.bed", sep="\t", index=False)

variant_dataframe #5,135

# filter for gnomad and get the population to be all
variant_dataframe_gnomad = variant_dataframe.query("study == 'GNOMAD_GENOMES'") #4,686
variant_dataframe_gnomad_ALL = variant_dataframe_gnomad.query("population == 'ALL'") #605
variant_dataframe_gnomad_ALL = variant_dataframe_gnomad_ALL.drop(['study', 'population'], axis=1) 
variant_dataframe_gnomad_ALL['INFO'] = 'AF=' + variant_dataframe_gnomad_ALL['INFO'].astype(str)
variant_dataframe_gnomad_ALL.to_csv("~/Documents/Projects/SexCheck/testing_data/myeloid/probes_myeloid_variant_gnomdALL.bed",
                                    sep="\t", index=False, )

# get any study but population=ALL
variant_dataframe_ALL = variant_dataframe.query("population == 'ALL'") 
variant_dataframe_ALL = variant_dataframe_ALL.drop(['study', 'population'], axis=1) 
variant_dataframe_ALL['INFO'] = 'AF=' + variant_dataframe_ALL['INFO'].astype(str)
variant_dataframe_ALL.to_csv("~/Documents/Projects/SexCheck/testing_data/myeloid/probes_myeloid_variant_ALL.bed",
                                    sep="\t", index=False, )

#via pycellbase
# variant_client does not mirror what is available on pycellbase :( 

# custom_config = {'rest': {'hosts': [host_address]}, 'version': 'v5', 'species': 'hsapiens'}
# cc = ConfigClient(custom_config)
# cbc = CellBaseClient(cc)
# cbc.show_configuration()['version']
# gv = cbc.get_variant_client()
# gv.get_help()