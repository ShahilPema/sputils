#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
import pandas as pd
import os

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infiles', type=click.Path(exists=True), required=True, nargs=-1)
@click.option('-o', '--outfile', type=click.Path(), required=True, help='Outfile.')
def main(infiles, outfile):
	"""
Concat CSV files.

\b
Example:
  pdconcatcsvs -o outfile.csv infile1.csv infile2.csv
	"""
	print('infiles: ', infiles)
	print('outfile: ', outfile)
	ins=[pd.read_csv(f) for f in infiles]
	tmp=pd.concat(ins)
	# (TODO) fix '.' without path
	os.makedirs(os.path.dirname(outfile), exist_ok=True)
	tmp.to_csv(outfile, index=False)

if __name__ == "__main__":
	main()
