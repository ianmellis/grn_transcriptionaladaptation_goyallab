import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import sem
import random

random.seed(10)
BOOTSTRAP_REPS = 1000

def filter_data(df, pos_only=False, sig_only=False):
    if 'Base_Mean' in df: # filter differently if DESeq Dataset
        # filter by base mean
        df = df[df['Base_Mean'] > 2.5]
        # add in support for p-value coloring
        df['significant'] = df['padj'].apply(lambda x: 1 if x<=0.05 else 0)
        # add in support for percent upreg
        df['upreg'] = df[['significant', 'FC2']].apply(lambda x: 1 if (x['FC2'] >= 0 and x['significant'] == 1) else 0, axis=1)
        # plot only KOs with significant paralog FC2s
        sig_only = True

    else: # standard rnaseq dataset
        # average points
        df = df.groupby(['KO_Gene', 'Paralog']).mean().reset_index()
        # filter by # raw read
        df = df[(df['Paralog_expr_count'] > 10) | (df['CTRL_expr_count'] > 10)]
        # filter by expression
        df = df[(df['Paralog_expr'] > 5) | (df['CTRL_expr'] > 5)]
        # add in support for percent upreg
        df['upreg'] = df['FC2'].apply(lambda x: 1 if x >= 0.5 else 0)

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

def gen_pvals(df, boot_dists):
    summarized = df.groupby('KO_Gene').mean().reset_index()
    data = []
    for _, row in summarized.iterrows():
        dist = boot_dists[row['KO_Gene']]
        perc_upreg_obs = row['upreg']
        perc_greater_than_obs = len(
            [i for i in dist if i >= perc_upreg_obs]) / len(dist)
        data.append([row['KO_Gene'], perc_greater_than_obs])

    df = pd.DataFrame(data, columns=['KO_Gene', 'p-value'])
        
    return df


def plot_data(df, boot_df, boot_dists, geo_id):
    fig_size = (13,8.5) 
    if 'Base_Mean' in df: # add in significance coloring
        hue = 'significant'
        palette = 'gist_gray_r'
        color = None
    else:
        color='black'
        hue = None
        palette = None
    
    df = filter_data(df, pos_only=False, sig_only=False)
    df = df.sort_values(by="KO_Gene")
    if len(set(df['KO_Gene'])) > 17:
        rotation = 45
    else:
        rotation = 0

    plt.figure(figsize=fig_size)
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                      font='sans-serif', font_scale=1, color_codes=True, rc=None)
    ax1 = sns.stripplot(x='KO_Gene', y='FC2', data=df, color=color, palette=palette, hue=hue)
    sns.barplot(x="KO_Gene", y="FC2", errorbar=None,color='gray', data=df, ax=ax1, alpha=.3)
    plt.xticks(rotation=rotation)
    plt.savefig('./analysis/FC2/' + geo_id + '.png', dpi=200)
    plt.clf()
    df.to_csv('./analysis/FC2/tabular_format/'+geo_id+'.csv')


    # setup df for percent upregulated
    pval_df = gen_pvals(df, boot_dists)
    boot_df = boot_df.groupby(['KO_Gene']).mean().reset_index()
    combined_df = boot_df.merge(df.groupby(['KO_Gene'])['upreg'].mean(), left_on='KO_Gene', right_on='KO_Gene')
    combined_df = combined_df.rename({'upreg': 'Paralogs', 'boot_upreg': 'Bootstrapped Genes'}, axis=1)
    upreg_data = pd.melt(combined_df, id_vars=['KO_Gene'], value_vars=['Paralogs', 'Bootstrapped Genes',])
    upreg_data = upreg_data.rename({'value':'Percent Upregulated'}, axis=1)
    upreg_data = upreg_data.merge(boot_df, how='left', left_on=['KO_Gene', 'Percent Upregulated'], right_on=['KO_Gene', 'boot_upreg'])
    upreg_data = upreg_data.merge(pval_df, how='left', left_on='KO_Gene', right_on='KO_Gene')

    # plot percent upregulated
    plt.figure(figsize=fig_size)
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                      font='sans-serif', font_scale=1, color_codes=True, rc=None)
    ax2 = sns.barplot(x="KO_Gene", y="Percent Upregulated", palette='Paired',
                      hue='variable', data=upreg_data)
    x_coords = [p.get_x() + 0.5*p.get_width() for p in ax2.patches]
    y_coords = [p.get_height() for p in ax2.patches]
    plt.errorbar(x=x_coords, y=y_coords, yerr=upreg_data["SE"], fmt="none", c="k")
    plt.xticks(rotation=rotation)
    plt.savefig('./analysis/percent_upregulated/' + geo_id + '.png', dpi=200)
    plt.clf()
    upreg_data.to_csv('./analysis/percent_upregulated/tabular_format/'+geo_id+'.csv')


    #plot knockout
    plt.figure(figsize=fig_size)
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                  font='sans-serif', font_scale=1, color_codes=True, rc=None)
    sns.stripplot(x='KO_Gene', y='KO_FC2', color='black',
                  jitter=False, data=df)
    plt.xticks(rotation=rotation)
    plt.savefig('./analysis/KO_FC2_plots/' + geo_id + '.png', dpi=200)
    plt.clf()


def process_deseq_data(deseq_file, paralog_data, ko_genes):
    #TODO: fix cases where there is a paralog with data, but padj is NaN
    fc_data = []
    deseq_file['baseMean'] = pd.to_numeric(deseq_file['baseMean'], errors='coerce')
    deseq_file.dropna(subset=['baseMean'], inplace=True)
    bootstrap_data = []
    bootstrap_dists = {}
    for gene in ko_genes:
        paralogs = list(paralog_data[paralog_data['Gene'] == gene]['Paralog'])
        perc_ids = list(paralog_data[paralog_data['Gene'] == gene]['Perc_ID'])

        gene_data = deseq_file[deseq_file['sampleKO']==gene]
        gene_data = gene_data.set_index('id')
        gene_data = gene_data.sort_values('baseMean')
        gene_data['rank'] = gene_data['baseMean'].rank()

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

        bootstrap_fc2_data = [[] for _ in range(len(paralogs))]
        bootstrap_padj_data = [[] for _ in range(len(paralogs))]
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

            rank = int(par_data['rank'][0])
            for _ in range(0,BOOTSTRAP_REPS):
                bootstrap_idx = random.randint(*random.choice([(rank-50, rank-1), (rank+1, rank+50)]))
                bootstrap_idx = 0 if (bootstrap_idx < 0) else bootstrap_idx
                bootstrap_idx = -1 if (bootstrap_idx >= len(gene_data)) else bootstrap_idx
                bootstrap_fc2_data[x].append(float(gene_data.iloc[bootstrap_idx]['log2FoldChange']))
                bootstrap_padj_data[x].append(float(gene_data.iloc[bootstrap_idx]['padj']))

            fc_data.append([gene, paralog, perc_ids[x], par_fc2, par_basemean, par_padj, ko_fc])
        
        if len(bootstrap_fc2_data) > 0:
            bootstrap_stats, bootstrap_dist = process_bootstrap_data(gene, gene, bootstrap_fc2_data, bootstrap_padj_data)
            bootstrap_data.append(bootstrap_stats)
            bootstrap_dists[gene] = bootstrap_dist

    df = pd.DataFrame(fc_data, columns=['KO_Gene', 'Paralog', 'Perc_ID', 'FC2', 'Base_Mean',
                      'padj', 'KO_FC2'])
    boot_df = pd.DataFrame(bootstrap_data, columns=['KO_Gene', 'Sample', 'boot_upreg', 'SE'])

    return df, boot_df, bootstrap_dists

def process_bootstrap_data(gene, sample, fc2_data, padj_data):
    # analyze the boostrapping samples
    # first remove paralogs that weren't found
    fc2_data = [l for l in fc2_data if len(l) != 0]
    padj_data = [l for l in padj_data if len(l) != 0]
    #if len(fc2_data) == 0:
    #    return [gene, sample, np.nan, np.nan]

    # start calculating averages
    perc_upregs = []
    for i in range(0, BOOTSTRAP_REPS):
        upreg_sum = 0
        for j, paralog_samples in enumerate(fc2_data):
            if len(padj_data) == 0: #TPM dataset
                if paralog_samples[i] > 0.5:
                    upreg_sum += 1
            else: #DESeq Data
                if paralog_samples[i] > 0 and padj_data[j][i] <= 0.05:
                    upreg_sum += 1

        perc_upreg = upreg_sum / len(fc2_data)
        perc_upregs.append(perc_upreg)

    avg_perc_upreg = sum(perc_upregs) / len(perc_upregs)
    se_perc_upreg = sem(perc_upregs)
    
    return [gene, sample, avg_perc_upreg, se_perc_upreg], perc_upregs 

def process_counts_data(norm_counts, counts, paralog_data, ko_genes, controls, condition):
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
    bootstrap_data = []
    bootstrap_dists = {}
    for gene in ko_genes:
        ko_ctrl_expr = pd.Series(norm_counts.loc[gene, 'control'])[0]

        relevant_samples = norm_counts.filter(
            regex='(^|[_])'+gene+'([_]|$)', axis=1)
        paralogs = list(paralog_data[paralog_data['Gene'] == gene]['Paralog'])
        perc_ids = list(paralog_data[paralog_data['Gene'] == gene]['Perc_ID'])

    
        # start calculating FC
        bootstrap_fc2_data = [[] for _ in range(len(paralogs))]
        for sample in relevant_samples:
            if condition != 'na' and condition not in sample:
                continue

            norm_counts = norm_counts.sort_values(sample)
            norm_counts['rank'] = norm_counts[sample].rank()
            
            # confirm KO
            ko_expr = pd.Series(norm_counts.loc[gene, sample])[0]
            ko_fc = np.log2(ko_expr/ko_ctrl_expr)

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
                
                paralog_expr = pd.Series(norm_counts.loc[paralog, sample])[0]
                paralog_expr_count = pd.Series(counts.loc[paralog, sample])[0]
    
                par_fc2 = np.log2(paralog_expr/par_ctrl_expr)

                rank = int(pd.Series(norm_counts.loc[paralog, 'rank'])[0])
                for _ in range(0,BOOTSTRAP_REPS):
                    bootstrap_idx = random.randint(*random.choice([(rank-50, rank-1), (rank+1, rank+50)]))
                    bootstrap_idx = 0 if (bootstrap_idx < 0) else bootstrap_idx
                    bootstrap_idx = -1 if (bootstrap_idx >= len(norm_counts)) else bootstrap_idx
                    sample_fc2 = np.log2(norm_counts.iloc[bootstrap_idx][sample] / norm_counts.iloc[bootstrap_idx]['control'])
                    bootstrap_fc2_data[x].append(sample_fc2)


                fc_data.append([gene, paralog, perc_ids[x],
                                sample, paralog_expr, paralog_expr_count, par_ctrl_expr, par_ctrl_expr_count, par_fc2, ko_fc])
            if len(bootstrap_fc2_data) > 0:
                bootstrap_stats, bootstrap_dist = process_bootstrap_data(gene, sample, bootstrap_fc2_data, [])
                bootstrap_data.append(bootstrap_stats)
                bootstrap_dists[gene] = bootstrap_dist

    df = pd.DataFrame(fc_data, columns=['KO_Gene', 'Paralog', 'Perc_ID', 'Sample', 'Paralog_expr', 'Paralog_expr_count', 'CTRL_expr', 'CTRL_expr_count','FC2', 'KO_FC2'])
    boot_df = pd.DataFrame(bootstrap_data, columns=['KO_Gene', 'Sample', 'boot_upreg', 'SE'])
    return df, boot_df, bootstrap_dists


def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')
    paralog_directory = './paralog_data/'
    norm_dataset_directory = './normalized_datasets/'
    raw_dataset_directory = './raw_datasets/'

    for x, row in dataset_metadata.iterrows():
        if x in [0,1,2]:
            continue

        condition = row['condition']
        deseq = row['DESeq']
        paralog_data = pd.read_csv(paralog_directory+row['GEO_ID']+'-paralogs.csv')
        ko_genes = row['ko_genes'].split(',')
        ctrl_samples = row['control_samples'].split(',')

        if deseq:
            seq_directory = norm_dataset_directory + 'DESeq/' + row['GEO_ID'] + '/' + 'differentialExpression_DESeq_allTargets.txt'
            seq_file = pd.read_csv(seq_directory, sep='\t')

            df, boot_df, boot_dists = process_deseq_data(seq_file, paralog_data, ko_genes)

        else:
            expression_data = pd.read_csv(norm_dataset_directory+row['GEO_ID']+'.csv')
            raw_expression_data = pd.read_csv(raw_dataset_directory+row['GEO_ID']+'.csv')
            
            df, boot_df, boot_dists = process_counts_data(expression_data, raw_expression_data, paralog_data, ko_genes, ctrl_samples, condition)

        plot_data(df, boot_df, boot_dists, row['GEO_ID'])


main()
