#!/usr/bin/env python3
import os
import re
import argparse
parser = argparse.ArgumentParser(
                    prog='Replaces any variable, x, written as ${x.inject} with new value y',
                    description='Inputs takens as x=y')
parser.add_argument('filein', type=str, help='Name fo the file to use as an input')
parser.add_argument('-f', '--force', action='store_true', help='Use flag to overwrite the output file')
parser.add_argument('-l', '--list', action='store_true', help='List the available arguments in the input file')
parser.add_argument('--sub', type=str, nargs='+', help='Space separated list in the format: VARIABLE_NAME=NEW_VALUE')
parser.add_argument('--output', type=str, help='Name of output file, required when using --sub')
args = parser.parse_args()

#args = argparse.Namespace(filein='../../../data/raw/plink.slurm.inject', list=True, sub=['JOB_NAME=my_new_name'], force=False, output='test.txt')

with open(args.filein, 'r') as f:
    text = f.read()

if args.list == True:
    all_available = re.findall('\$\{.*.inject\}', text)
    all_available = [x.removeprefix('${').removesuffix('.inject}') for x in all_available]
    all_available = list(set(all_available)) # Remove duplicates
    print(f"Variables to replace in {args.filein} :")
    print(all_available)

if args.sub is not None:
    # Don't continue if no output is given
    if args.output is None:
        raise Exception('No output file given! Required when using --sub')
        
    # Replace variables if appropriate
    to_replace = {f"${{{x.split('=')[0]}.inject}}":x.split('=')[1] for x in args.sub}
    for key, item in to_replace.items():
        text = text.replace(key, item)
        
    # Write out the file
    if os.path.isfile(args.output) and args.force==False:
        raise Exception('Output file already exists, use -f to force rewrite')
    with open(args.output, 'w') as f:
        f.write(text)