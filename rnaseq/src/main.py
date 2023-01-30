import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def filter_data(df, pos_only=False, sig_only=False):
    if 'Base_Mean' in df: # filter differently if DESeq Dataset
        # filter by base mean
        df = df[df['Base_Mean'] > 2.5]
        # add in support for p-value coloring
        df['significant'] = df['padj'].apply(lambda x: 1 if x<=0.05 else 0)

    else: # standard rnaseq dataset
        # average points
        df = df.groupby(['KO_Gene', 'Paralog']).mean().reset_index()
        # filter by # raw read
        df = df[(df['Paralog_expr_count'] > 10) | (df['CTRL_expr_count'] > 10)]
        # filter by expression
        df = df[(df['Paralog_expr'] > 5) | (df['CTRL_expr'] > 5)]

    if pos_only:
        # filter out any average negative expression
        summarized = df.groupby('KO_Gene').mean()
        relevant_samples = summarized[summarized['FC2'] > 0.1]
        df = df[df['KO_Gene'].isin(list(relevant_samples.index))]

    if sig_only:
        # filter out any average non-significant expression
        summarized = df.groupby('KO_Gene').sum()
        relevant_samples = summarized[summarized['significant'] > 0]
        df = df[df['KO_Gene'].isin(list(relevant_samples.index))]

    return df

def plot_data(df, condition):
    #TODO: plot fraction of paralogs upreg: like other paper
    if condition != 'na':
        df = df[df['Sample'].str.contains(condition)]
    
    df = filter_data(df, pos_only=False, sig_only=True)

    df = df.sort_values(by="KO_Gene")
    #fig, ax1 = plt.subplots()
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                      font='sans-serif', font_scale=1, color_codes=True, rc=None)

    
    ax1 = sns.stripplot(x='KO_Gene', y='FC2', data=df, palette='gist_gray_r', hue='significant')
    sns.barplot(x="KO_Gene", y="FC2", errorbar=None,color='gray', data=df, ax=ax1, alpha=.3)
    plt.show()

    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                  font='sans-serif', font_scale=1, color_codes=True, rc=None)

    sns.stripplot(x='KO_Gene', y='KO_FC2', color='black', jitter=False, data=df)
    plt.show()

    #sns.set_theme(context='notebook', style='whitegrid', palette='deep',
    #              font='sans-serif', font_scale=1, color_codes=True, rc=None)

    #df.plot(kind='bar', stacked=True, x='KO_Gene', y='KO_FC2')
    #plt.show()


def process_deseq_data(deseq_file, paralog_data, ko_genes):
    fc_data = []

    for gene in ko_genes:

        paralogs = list(paralog_data[paralog_data['Gene'] == gene]['Paralog'])
        perc_ids = list(paralog_data[paralog_data['Gene'] == gene]['Perc_ID'])

        gene_data = deseq_file[deseq_file['sampleKO']==gene]
        gene_data = gene_data.set_index('id')

        # confirm KO
        ko_data = gene_data.filter(regex='(^|[_])'+gene+'([_]|$)', axis=0)
        
        if len(ko_data) > 1:
            ko_data = ko_data.replace("NA", np.nan)
            ko_data = ko_data.dropna()
        if len(ko_data) != 1:
            ko_fc = np.nan
            print('Warning: ' + gene + ' KO gene not found')
        else:
            ko_fc = float(ko_data['log2FoldChange'][0])


        for x, paralog in enumerate(paralogs):
            try:
                par_data = gene_data.filter(regex='(^|[_])'+paralog+'([_]|$)', axis=0)
                if len(par_data) > 1:
                    par_data = par_data.replace("NA", np.nan)
                    par_data = par_data.dropna()
                assert len(par_data) == 1
            except:
                print('Warning: ' + paralog + ' paralog not found')
                continue

            par_fc2 = float(par_data['log2FoldChange'][0])
            par_basemean = float(par_data['baseMean'][0])
            par_padj = float(par_data['padj'][0])

            fc_data.append([gene, paralog, perc_ids[x], par_fc2, par_basemean, par_padj, ko_fc])

    df = pd.DataFrame(fc_data, columns=['KO_Gene', 'Paralog', 'Perc_ID', 'FC2', 'Base_Mean', 'padj','KO_FC2'])

    return df


def process_counts_data(norm_counts, counts, paralog_data, ko_genes, controls):
    norm_counts = norm_counts.set_index('gene_name')
    counts = counts.set_index('gene_name')

    # calulate control expression
    control_samples = [list(norm_counts.filter(regex='^'+x+'$', axis=1)) for x in controls]
    control_samples = list(itertools.chain.from_iterable(control_samples))
    norm_counts['control'] = norm_counts[control_samples].mean(axis=1)
    counts['control'] = counts[control_samples].mean(axis=1)

    norm_counts += 1

    
    #TODO sum the expression when there are multiple genes found
    fc_data = []
    for gene in ko_genes:
        ko_ctrl_expr = pd.Series(norm_counts.loc[gene, 'control'])[0]

        relevant_samples = norm_counts.filter(
            regex='(^|[_])'+gene+'([_]|$)', axis=1)
        paralogs = list(paralog_data[paralog_data['Gene'] == gene]['Paralog'])
        perc_ids = list(paralog_data[paralog_data['Gene'] == gene]['Perc_ID'])

        for x, paralog in enumerate(paralogs):
            # get control expression:
            try:
                par_ctrl_expr = pd.Series(
                    norm_counts.loc[paralog, 'control'])[0]
                par_ctrl_expr_count = pd.Series(
                    counts.loc[paralog, 'control'])[0]
            except KeyError:
                print('Warning: ' + paralog + ' not found')
                continue
            
            # start calculating FC
            for sample in relevant_samples:
                # confirm KO
                ko_expr = pd.Series(norm_counts.loc[gene, sample])[0]
                ko_fc = np.log2(ko_expr/ko_ctrl_expr)

                paralog_expr = pd.Series(norm_counts.loc[paralog, sample])[0]
                paralog_expr_count = pd.Series(counts.loc[paralog, sample])[0]

                par_fc2 = np.log2(paralog_expr/par_ctrl_expr)
                fc_data.append([gene, paralog, perc_ids[x],
                               sample, paralog_expr, paralog_expr_count, par_ctrl_expr, par_ctrl_expr_count, par_fc2, ko_fc])

    df = pd.DataFrame(fc_data, columns=['KO_Gene', 'Paralog', 'Perc_ID', 'Sample', 'Paralog_expr', 'Paralog_expr_count', 'CTRL_expr', 'CTRL_expr_count','FC2', 'KO_FC2'])

    return df


def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')
    paralog_directory = './paralog_data/'
    norm_dataset_directory = './normalized_datasets/'
    raw_dataset_directory = './raw_datasets/'

    for x, row in dataset_metadata.iterrows():
        if x in [0,2]:
            continue

        condition = row['condition']
        deseq = row['DESeq']
        paralog_data = pd.read_csv(paralog_directory+row['GEO_ID']+'-paralogs.csv')
        ko_genes = row['ko_genes'].split(',')
        ctrl_samples = row['control_samples'].split(',')

        if deseq:
            seq_directory = norm_dataset_directory + 'DESeq/' + row['GEO_ID'] + '/' + 'differentialExpression_DESeq_allTargets.txt'
            seq_file = pd.read_csv(seq_directory, sep='\t')

            df = process_deseq_data(seq_file, paralog_data, ko_genes)

        else:
            expression_data = pd.read_csv(norm_dataset_directory+row['GEO_ID']+'.csv')
            raw_expression_data = pd.read_csv(raw_dataset_directory+row['GEO_ID']+'.csv')
            
            df = process_counts_data(expression_data, raw_expression_data, paralog_data, ko_genes, ctrl_samples)

        plot_data(df, condition)


main()
