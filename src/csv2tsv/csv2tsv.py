#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
import pandas as pd
import os
import csv

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infile', type=click.Path(exists=True), required=True, nargs=1)
@click.option('-o', '--outfile', type=click.Path(), required=True, help='Outfile.')
def main(infile, outfile):
	"""
Convert csv to tsv files.

\b
Example:
  csv2tsv -o outfile.csv  infile2.csv
	"""
	print(f"Converting {os.path.basename(infile)} to tsv")
	os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
	with open(infile, 'r') as csv_in, open(outfile, 'w') as tsv:
		csv_in = csv.reader(csv_in, delimiter=',', skipinitialspace=True)
		tsv = csv.writer(tsv, delimiter="\t")
		for row in csv_in:
			tsv.writerow(row)

if __name__ == "__main__":
	main()
