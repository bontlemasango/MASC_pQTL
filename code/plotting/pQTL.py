import os
import re
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

parser = argparse.ArgumentParser(
                    prog='Script to make an individual page of manahattan and qq plots for proteins',
                    description='')
parser.add_argument('folder', type=str, help='Folder containing GWAS file outputs. Files should be int he format chrX_protein_name.regenie')
parser.add_argument('hgnc_file', type=str, help='hgnc tsv containing the gene locations')
parser.add_argument('--genename', action='store_true', help='Use if protein IDs are given as gene names instead of uniprot IDs')
parser.add_argument('--output', type=str, help='Output file name', default='pQTL_plot.png')
args = parser.parse_args()#['/data/gen1/Brandon_PhD/proteomic_GWAS_freeze2/data/processed/regenie_step2/', '/data/gen1/Brandon_PhD/proteomic_GWAS_freeze2/code/plotting/tmp_data/hgnc_ucsc_b38_genes.tsv'])


all_files = os.listdir(args.folder)
all_proteins = pd.DataFrame(dict(filename=all_files))
all_proteins = (all_proteins
        .query('filename!="log_files"'))
all_proteins = (all_proteins
        .query('~filename.str.endswith(".log")')
                .assign(chromosome=all_proteins['filename'].str.split('_').str[1].replace('X', 23).astype(int))
                .assign(protein_id=all_proteins['filename'].str.replace('^chr_.*?_', '', regex=True).replace('.regenie$','', regex=True))
                .sort_values(['protein_id','chromosome'])
               )
all_proteins['chromosome'] = all_proteins['chromosome'].astype(int)

hgnc = pd.read_csv(args.hgnc_file, sep='\t', usecols=['chrom', 'chromStart','symbol','geneName','uniprot_ids','mane_select'])
hgnc = hgnc.query('~uniprot_ids.isna()')
hgnc ['symbol'] = hgnc ['symbol'].str.replace('-','.')

mydir = lambda x: os.path.join(args.folder, x)
def read_dfs(protein_id):
    protein = all_proteins.query('protein_id==@protein_id')
    dfs = [pd.read_csv(mydir(x), sep=' ', usecols=['CHROM','GENPOS', 'LOG10P']).query('LOG10P>3') for x in protein['filename']]
    df = pd.concat(dfs)
    df['PROTEIN'] = protein_id
    return df


def merge_datasets(gwas_df, hgnc_df):
    merge_on = 'symbol' if args.genename else 'uniprot_ids'
    hgnc_df = hgnc_df.copy()
    hgnc_df['cognate_chrom'] = hgnc_df['chrom'].str.removeprefix('chr').replace('X',23)
    hgnc_df['cognate_pos'] = hgnc_df['chromStart']
    hgnc_df = hgnc_df[['cognate_chrom', 'cognate_pos', 'symbol', 'uniprot_ids']]
    out = gwas_df.copy()
    out['PROTEIN'] = out['PROTEIN'].str.split('_').str[0]
    out = out.merge(hgnc_df, how='left', left_on='PROTEIN', right_on=merge_on)
    out = out.rename(columns={'CHROM':'chrom','GENPOS':'pos',
                              'PROTEIN':'protein_group','symbol':'gene',
                              'LOG10P':'log10p',
                             'cognate_chrom':'chrom_gene', 'cognate_pos':'transcript_start'})
    out['chrom_gene'] = pd.to_numeric(out['chrom_gene'], errors='coerce')
    return out


def calculate_x_axis(df):
    df = df.copy()
    minpos_dict = df.groupby('chrom')['pos'].min().to_dict()
    df['relative_pos'] = df['pos'] - df['chrom'].map(minpos_dict)
    
    maxpos_dict = df.groupby('chrom')['relative_pos'].max().cumsum().shift(1).fillna(0).to_dict()
    df['x_axis'] = df['relative_pos'] + df['chrom'].map(maxpos_dict)
    df['y_axis'] = df['transcript_start'] - df['chrom_gene'].map(minpos_dict) + df['chrom_gene'].map(maxpos_dict)
    df = df.drop(columns=['relative_pos'])
    #df['y_axis'] = df['chrom_gene'].map(maxpos_dict) - df['chrom_gene'].map(minpos_dict)
    #df['y_axis'] = df['transcript_start'] + df['y_axis']
    return df


# ## Get a list of the lead SNPs

def lead_snps(df, max_dist=1_000_000):
    df = df.copy()
    df = df.sort_values('log10p', ascending=False).reset_index(drop=True)
    for i in range(df.shape[0]):
        chrom_match = df['chrom']==df.loc[i, 'chrom']
        dist = abs(df['pos']-df.loc[i, 'pos'])
        df = df.loc[~(chrom_match&(dist<max_dist)&(df.index!=i)), :].reset_index(drop=True)
        if i>=df.shape[0]-1:
            break
    return df


# ## Plotting

def create_base(df):
    fig, axs = plt.subplots(2, 1, height_ratios=(1,2), figsize=(20,18))
    fig.subplots_adjust(hspace=0.05)
    
    # Set the appropriate limits and tick positions for the x-axis of both plots
    xlims = [df['x_axis'].min(), df['x_axis'].max()]
    xticks = df.groupby('chrom')['x_axis'].mean()
    for ax in axs:
        ax.set_xlim(xlims)
        ax.set_xticks(xticks, minor=False)
        ax.set_xticklabels(xticks.index, minor=False)
    
    # Set y-axis for manhattan
    axs[0].set_ylim([df['log10p'].min(), df['log10p'].max()+1])
    
    # Set y-axis for pQTL plot
    axs[1].set_ylim(xlims)
    axs[1].set_yticks(xticks, minor=False)
    axs[1].set_yticklabels(xticks.index, minor=False)
    
    # Add gridlines for each chromosome
    axs[1].set_xticks(df.groupby('chrom')['x_axis'].max(), minor=True)
    axs[1].xaxis.grid(True, which='minor')
    axs[1].set_yticks(df.groupby('chrom')['x_axis'].max(), minor=True)
    axs[1].yaxis.grid(True, which='minor')
    
    return fig, axs


def rockefeller(df, ax, sub_pthresh=5*10**-6, pthresh=5*10**-8,
                text_thresh=None, highlight_leadsnps=False,
                markersize=8, leadsnp_dist=1_000_000,
               fontsize=12):
    
    plot_colours = df['chrom'].map(lambda x: 'black' if x%2==1 else 'grey')
    # Calculate the lead SNPs and everythin within 1Mb
    if highlight_leadsnps:
        lead_snps_df = lead_snps(df[df['log10p']>-np.log10(pthresh)], max_dist=leadsnp_dist)
        for i in range(lead_snps_df.shape[0]):
            plot_colours[(abs(df['x_axis']-lead_snps_df.loc[i, 'x_axis'])<=leadsnp_dist)&(df['chrom']==lead_snps_df.loc[i, 'chrom'])] = 'green'
            
    # Plot non-green points first
    myfilt = plot_colours!='green'
    ax.scatter(df.loc[myfilt, 'x_axis'], df.loc[myfilt, 'log10p'], c=plot_colours[myfilt], s=markersize)
    ax.scatter(df.loc[~myfilt, 'x_axis'], df.loc[~myfilt, 'log10p'], c=plot_colours[~myfilt], s=markersize)
    # Add significance lines
    ax.hlines([-np.log10(sub_pthresh), -np.log10(pthresh)],
              xmin=plot_df['x_axis'].min(), xmax=plot_df['x_axis'].max(),
              color=['green', 'red'], linestyles=['--','-'], alpha=0.5)
    
    # Highlight the lead snp gene names
    if text_thresh is not None:
        text_df = lead_snps_df.loc[lead_snps_df['log10p']>-np.log10(text_thresh), ['x_axis', 'log10p', 'gene']] 
        plot_text = [ax.text(d['x_axis'], d['log10p'], d['gene'], fontsize=fontsize) for _,d in text_df.iterrows()]
        #adjustText.adjust_text(plot_text, ax=ax, force_text=(0.1,0.2), force_static=(0,3), time_lim=2, arrowprops=dict(arrowstyle='-'))


def pQTL(df, ax, sub_pthresh=5*10**-6, pthresh=5*10**-8,
         guidelines=False, markersize=8, leadsnp_dist=1_000_000):
    
    # Filter the dataframe to only take relevant p-values
    lead_snps_df = lead_snps(df[df['log10p']>-np.log10(sub_pthresh)], max_dist=leadsnp_dist)
    
    # Plot the sub-genome wide significant pQTLs
    lead_snps_df['colour'] = lead_snps_df['chrom'].map(lambda x: 'grey' if x%2==1 else 'silver')
    lead_snps_df.loc[lead_snps_df['chrom']==lead_snps_df['chrom_gene'], 'colour'] = 'green'
    
    # Modify markers for above the pthresh and below pthresh
    lead_snps_df['marker'] = lead_snps_df['log10p'].map(lambda x: 'o' if x>-np.log10(pthresh) else 'x')
    
    for marker, d in lead_snps_df.groupby('marker'):
        ax.scatter(d['x_axis'], d['y_axis'], c=d['colour'], marker=marker, alpha=0.6)
        
    # Add guidelines if necessary
    if guidelines:
        # Add the x=y line
        ax.axline([df['x_axis'].min()]*2, slope=1, color='red', linewidth=1, alpha=0.2)
        # Add in vertical lines to make spotting the association easier
        vlines = lead_snps_df.loc[lead_snps_df['marker']=='o', ['x_axis', 'y_axis']]
        vlines['y_max'] = df['y_axis'].max()
        ax.vlines(vlines['x_axis'], ymin=vlines['y_axis'], ymax=vlines['y_max'], color='green', linestyles=':', alpha=0.3)


gwas_df = []
for i, protein in enumerate(all_proteins['protein_id'].unique()):
    gwas_df.append(read_dfs(protein))
    if i>1:
        break
gwas_df = pd.concat(gwas_df)

plot_df = merge_datasets(gwas_df, hgnc).pipe(calculate_x_axis)

# +
gwas_thresh = 5*10**-8
sub_gwas_thresh = 5*10**-6

myfig, myaxs = create_base(plot_df)
rockefeller(plot_df, myaxs[0],
            sub_pthresh=sub_gwas_thresh, pthresh=gwas_thresh,
            text_thresh=gwas_thresh/100, highlight_leadsnps=True, markersize=3)
pQTL(plot_df, myaxs[1], 
     sub_pthresh=sub_gwas_thresh, pthresh=gwas_thresh,
     markersize=3, guidelines=True)
# -
myfig.savefig(args.output)

