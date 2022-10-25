#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
import pandas as pd
import os

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infiles', type=click.Path(exists=True), required=True, nargs=-1)
@click.option('--ids_header', '-n', type=click.STRING, default='File ID', required=False)
@click.option('--ids', '-i', multiple=True, required=False)
@click.option('-o', '--outfile', type=click.Path(), required=True, help='path/to/outfile.csv')
def main(infiles, ids, outfile, ids_header):
    """
Concat CSV files.

\b
Example:
  pdconcatcsvs -o path/to/outfile.csv infile1.csv infile2.csv
  OR
  pdconcatcsvs -o path/to/outfile.csv -i id1 -i id2 infile1.csv infile2.csv
  OR
  pdconcatcsvs -o path/to/outfile.csv -n 'Sample IDs' -i id1 -i id2 infile1.csv infile2.csv
    """
    ins=[pd.read_csv(f) for f in infiles]
    if len(ids)>0:
        if len(infiles) != len(ids):
            raise Exception("The number of sample IDs must match the number of infiles")

        in_file_df_list=[]
        for id, df in zip(ids, ins):
            df.insert(0, ids_header, id)
            in_file_df_list += [df]
        ins=in_file_df_list

    all_smpls_df=pd.concat(ins, ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    all_smpls_df.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
