import os
import re
import argparse
from os import path

# +
# Variables to reference
SUFFIX = '.regenie'
SEP = ' '
LOG10P_COLUMN = 'LOG10P'

# Used for the commandline inputs
parser = argparse.ArgumentParser(
    prog='SNP filter - log10P',
    description='Extract only the SNPs with -log10P >= the minimum value given',
    usage=f'''
    Example:
    myfolder/results -m 4.5 -o all_results.tsv
    Will go through every file in myfolder/results that ends with {SUFFIX},
    read it using pandas (assuming tab-delimited),
    add a column with the protein name,
    filter the results for every SNP with -log10P more than (or equal to) the value given for -m,
    and write it out to all_results.tsv
    '''
)
parser.add_argument('folder', help=f'folder containing {SUFFIX} files to filter')
parser.add_argument('-m', '--min-log10p', type=float, default=7.3, help='Minimum -log10P to filter results by. Deafults to genome wide significance threshold: 5E-8')
parser.add_argument('-o', '--output', default='MASC_pQTLs.{min_log10p}.tsv', help='Output file name.')
#args = parser.parse_args('../tst/chrom_split -m 4'.split())
args = parser.parse_args()

# Used to record the min-log10p input into the program
args.output = args.output.replace('{min_log10p}', str(args.min_log10p).replace('.','_'))


# -

def list_files(folderpath):
    ''' List all the files in the folder that end with .regenie '''
    out = [path.join(folderpath, x) for x in os.listdir(folderpath) if x.endswith(SUFFIX)]
    return out


def get_proteinname(filename):
    ''' Figure out the protein name based on the filename ie: chr_2_[PROTEIN_NAME].regenie'''
    protein_name = path.basename(filename)
    protein_name = protein_name.removesuffix(SUFFIX)
    protein_name = re.sub('^chr_[0-9]+_', '', protein_name)
    return protein_name


def filter_file(filename, add_protein_column=True):
    ''' Get all the lines that have LOG10P>min_log10p (from commandline arguments)'''
    out = []
    with open(filename, 'r') as f:
        header = next(f).strip().split(SEP)
        for line in f:
            data = line.strip().split(SEP)
            if float(data[header.index(LOG10P_COLUMN)]) > args.min_log10p:
                out.append(line)
    if add_protein_column:
        out = [x.strip('\n') + SEP + get_proteinname(filename) for x in out]
    return out


def update_header(filename):
    with open(filename, 'r') as f:
        header = next(f).strip('\n')
        header = f'{header}{SEP}PROTEIN_ID\n'
    return header


if __name__=='__main__':
    with open(args.output, 'w') as f_out:
        # List the files
        all_files = list_files(args.folder)
        # Add a header first
        f_out.write(update_header(all_files[0]))
        # For each file, write the results out to your output file
        for i, file in enumerate(all_files):
            print(f"{(i+1)/len(all_files)*100: .2f}", end='\r')
            results = filter_file(file) + ['']
            f_out.write('\n'.join(results))
        print('Done!')



