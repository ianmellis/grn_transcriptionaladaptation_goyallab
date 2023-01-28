import pandas as pd
import requests
import sys
from gprofiler import GProfiler


def get_paralog_info(gene_name, species) -> list:
    server = "https://rest.ensembl.org"
    ext = "/homology/symbol/"+species+"/" + gene_name + "?"

    r = requests.get(server+ext, headers={"Content-Type": "application/json",
                     "format": "condensed", "type": "paralogues", 'target_species': species})

    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()

    # parse data into list of paralog ENSEMBL IDs
    paralogs = []
    similarities = []
    for matched_gene in decoded['data']:
        for paralog_data in matched_gene['homologies']:
            ens_id = paralog_data['target']['id']
            if ('ENSG0' in ens_id and species=='human') or ('ENSMUSG0' in ens_id and species=='mouse'):
                paralogs.append(ens_id)
                similarities.append(paralog_data['target']['perc_id'])

    # converts IDs to gene names
    if len(paralogs) > 0:
        g_species = 'hsapiens' if species == 'human' else 'mmusculus'
        gp = GProfiler(return_dataframe=True)
        converted_paralogs = gp.convert(organism=g_species,
                                        query=paralogs,
                                        target_namespace='ENSG')
        cleaned_paralogs = [x for x in list(converted_paralogs['name']) if x!='None']
        return cleaned_paralogs, similarities
    else:
        return paralogs, similarities


def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')

    for index, row in dataset_metadata.iterrows():
        if index in [0]:
            continue
        species = row['species']
        ko_genes = row['ko_genes'].split(',')

        paralog_data = []
        for i, gene in enumerate(ko_genes):
            paralogs, similarities = get_paralog_info(gene, species)
            for x, paralog in enumerate(paralogs):
                paralog_data.append([gene, paralog, similarities[x]])
                print(gene, paralog)

            print(str(i+1)+' of ' + str(len(ko_genes)) + ' ko_genes complete')
        df = pd.DataFrame(paralog_data, columns=['Gene', 'Paralog', 'Perc_ID'])
        #df = pd.DataFrame.from_dict(paralog_data, orient='index')
        #df = df.transpose()
        df.to_csv('./paralog_data/'+row['GEO_ID']+'-paralogs.csv')

        print(str(index+1)+' of ' +
              str(len(dataset_metadata)) + ' datasets complete')


main()

#TODO set ENSEMBL Version
#TODO explore "None" in DEAF1