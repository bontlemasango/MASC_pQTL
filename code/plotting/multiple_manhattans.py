import os
import re
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')

parser = argparse.ArgumentParser(
                    prog='Script to make an individual page of manahattan and qq plots for proteins',
                    description='')
parser.add_argument('folder', type=str, help='Folder containing GWAS file outputs. Files should be int he format chrX_protein_name.regenie')
parser.add_argument('page_number', type=int, help='page number to create')
parser.add_argument('--output-prefix', type=str, help='prefix for the output file. Outputs have suffix .pageX.pdf', default='plot')
parser.add_argument('--max-pages', action='store_true', help='Prints the number of pages needed for everything (no output is made)')
parser.add_argument('--nrows', type=int, default=5, help='Number of rows per page')
args = parser.parse_args()#(['/data/gen1/Brandon_PhD/proteomic_GWAS_freeze2/data/processed/regenie_step2/', '3', '--output-prefix', 'tst_plot'])

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


if args.max_pages:
    print(f"Number of pages needed to print all manhattan plots: {np.ceil(all_proteins['protein_id'].nunique()/args.nrows): .0f}")
    exit()
else:
    NROWS = args.nrows # Number of rows per page
    PAGE = args.page_number
    OUTPUT = f"{args.output_prefix}.page{args.page_number}.pdf"
    proteins_to_take = NROWS*(PAGE-1)
    proteins_to_take = (all_proteins['protein_id']
                        .drop_duplicates()
                        .iloc[proteins_to_take:min(proteins_to_take+NROWS, all_proteins['protein_id'].nunique())]
                       )
    plot_proteins = all_proteins.query('protein_id.isin(@proteins_to_take)').reset_index(drop=True)

# +
MEDIAN_CHISQ = chi2.ppf(0.5, 1)

mydir = lambda x: os.path.join(args.folder, x)

def read_dfs(protein_id):
    protein = plot_proteins.query('protein_id==@protein_id')
    dfs = [pd.read_csv(mydir(x), sep=' ', usecols=['CHROM','GENPOS','CHISQ', 'LOG10P']) for x in protein['filename']]
    df = pd.concat(dfs)
    return df


# +
def calc_lambdagc(df):
    return df['CHISQ'].median()/MEDIAN_CHISQ

def plot_qq(df, ax=None):
    ''' Note: this uses the fact that he index is NOT reset ie multiple 0s'''
    lambdagc = calc_lambdagc(df)
    df = df.loc[:,['LOG10P']].copy().sort_values('LOG10P').reset_index(drop=True)
    df['expected'] = -np.log10(np.linspace(1, 0, df.shape[0], endpoint=False))
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(5,5))
    ax.axline((0, 0), slope=1, color='red', alpha=0.5)
    ax.scatter(df['expected'], df['LOG10P'], color='black', s=3, rasterized=True)
    ax.text(0.01, 0.99,  f"lambda_gc={lambdagc: .5f}", ha='left', va='top', fontsize=11, transform = ax.transAxes, usetex=False)

def plot_manhattan(df, ax=None, sig_line=-np.log10(5*10**-8/(all_proteins['protein_id'].nunique())), subsig_line=-np.log10(5*10**-8)):
    df = df.query('LOG10P>1.2').copy()
    df['CHROM'] = pd.to_numeric(df['CHROM'].replace({'X':23, 'Y':24}), errors='coerce').astype('Int64')
    df['x'] = df.groupby('CHROM')['GENPOS'].transform('min')
    df['x'] = df['GENPOS'] - df['x']
    df['x'] += 1_000_000
    df['x'] += df['CHROM'].map(df.groupby('CHROM')['x'].max().shift().cumsum().fillna(0).to_dict())
    df['c'] = (df['CHROM']%2).map({0:'grey',1:'black'})
    
    x_labels = df.groupby('CHROM')['x'].mean()
    x_labels_idx = x_labels.values.reshape(-1)
    x_labels = x_labels.index
    
    if ax is None:
        fix, ax = plt.subplots(figsize=(8, 5))
    
    ax.set_xticks(x_labels_idx)
    ax.set_xticklabels(x_labels)
    ax.scatter(df['x'], df['LOG10P'], c=df['c'], s=3, rasterized=True)
    ax.hlines([subsig_line, sig_line],
              xmin=df['x'].min(), xmax=df['x'].max(), alpha=0.5, colors=['green', 'red'], linestyles=['--', '-'])
    ax.set_ylim(df['LOG10P'].min(), max(df['LOG10P'].max(), sig_line)+1)
    ax.set_xlim(df['x'].min(), df['x'].max())


# +
# Create and write the file as necessary
if os.path.isfile(OUTPUT):
    raise Exception(f"{OUTPUT} already exists, delete to make new pdf")

pdf = matplotlib.backends.backend_pdf.PdfPages(OUTPUT)

myscale = 1
fig, axs = plt.subplots(ncols=2, nrows=NROWS, width_ratios=[2,1], figsize=(8.27*myscale, 11.69*myscale), rasterized=False)
fig.subplots_adjust(hspace=0.4)
for ax_idx, protein_id in enumerate(proteins_to_take.to_list()):
    plot_title = protein_id
    df = read_dfs(protein_id)
    axs[ax_idx, 0].set_title(plot_title, loc='left', y=1.0)
    plot_manhattan(df, axs[ax_idx, 0])
    plot_qq(df, axs[ax_idx, 1])
    print(f"{ax_idx+1}/{NROWS}", end='\n')
    
pdf.savefig(fig, dpi=150)
pdf.close()
print('Done!')
# -


