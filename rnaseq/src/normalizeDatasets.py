import itertools
import pandas as pd
import numpy as np
from bioinfokit.analys import norm
import re


def add_gene_lengths(df, species):

    if species == 'human':
        annot = './annotations/Homo_sapiens.GRCh38.108.gtf'
    if species == 'mouse':
        annot = './annotations/Mus_musculus.GRCm39.108.gtf'

    wantedGenes = set(df['gene_name'])
    length_dict = {}

    with open(annot, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split('\t')
                annotations = fields[8].split(' ')
                geneName = re.sub('[";]', '', annotations[5])
                if geneName in wantedGenes and fields[2] == 'gene':
                    length_dict[geneName] = int(
                        int(fields[4]) - int(fields[3]))

    df['length'] = df['gene_name'].map(length_dict)

    return df


def normalize(counts, species):

    # normalize raw counts using CPM method
    counts = add_gene_lengths(counts, species)
    counts = counts.set_index('gene_name')

    nm = norm()
    nm.tpm(df=counts, gl='length')
    norm_counts = nm.tpm_norm

    return norm_counts


def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')
    dataset_directory = './raw_datasets/'
    normalized_directory = './normalized_datasets/'

    for x, row in dataset_metadata.iterrows():
        expression_data = pd.read_csv(dataset_directory+row['GEO_ID']+'.csv')

        df = normalize(expression_data, row['species'])
        df.to_csv(normalized_directory+row['GEO_ID']+'.csv')
        print(str(x+1)+' of ' + str(len(dataset_metadata)) + ' datasets complete')


main()
