#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
import pandas as pd
import os
import re
import sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infiles', type=click.Path(exists=True), required=False, nargs=-1)
@click.option('-i', '--ids', type=click.STRING, required=False, help='Provide string of comma separated sample IDs.')
@click.option('-d', '--indir', type=click.Path(), required=False, help='path/to/dir containing infiles.csv')
@click.option('-o', '--outfile', type=click.Path(), required=True, help='path/to/outfile.csv')
def main(infiles, ids, indir, outfile):
    """
Concat CSV files.

\b
Example:
  pdconcatcsvs -o paht/to/outfile.csv -i 'id1,id2' infile1.csv infile2.csv
  OR
  pdconcatcsvs -o path/to/outfile.csv -d 10x_csvs/
    """
    in_file_df_list = []    
    in_subfolders = os.listdir(indir)
    if indir:
        for subfolder in in_subfolders:
            smpl_id = re.search(r'\w\d{3}', subfolder)[0] #match ex. D001 from subfolder name
            id_df = pd.DataFrame(data={'Sample ID':[smpl_id]})
            qc_sum = pd.read_csv(f"{indir}/{subfolder}/metrics_summary.csv")
            qc_sum_w_id = pd.concat([id_df, qc_sum], ignore_index=False, sort=False, axis=1)
            in_file_df_list += [qc_sum_w_id]     
    elif infiles:
        try:
            id_list = ids.split(',')
            id_list = [x.strip() for x in id_list]
            if len(infiles) == len(id_list):
                count = 0
                for file in infiles:
                    id_df = pd.DataFrame(data={'Sample Id':[id_list[count]]})
                    qc_sum = pd.read_csv(file)
                    qc_sum_w_id = pd.concat([id_df, qc_sum], ignore_index=False, sort=False, axis=1)
                    in_file_df_list += [qc_sum_w_id]
                    count += 1
            else:
                raise Exception("The number of sample IDs must match the number of infiles")
        except Exception as e:
            print("No sample IDs were provided.")
            sys.exit(-1)
    all_smpls_df=pd.concat(in_file_df_list, ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    all_smpls_df.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
