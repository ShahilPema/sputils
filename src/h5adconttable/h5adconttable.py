#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
import pandas as pd
import scanpy as sc
import numpy as np

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infile', type=click.Path(exists=True), required=True, nargs=1)
@click.option('-x', '--xvariable', type=click.STRING, required=True, help='Variable to be used as row index.')
@click.option('-y', '--yvariable', type=click.STRING, required=True, help='Variabel to be used as columns.')
@click.option('-o', '--outfile', type=click.Path(), required=False, help='Outfile.')
@click.option('-b', '--bname', type=click.STRING, required=False, help='Basename.')
def main(infile, xvariable, yvariable, outfile, bname):
	"""
Create contingency table from variable in an h5ad file.

\b
Example:
  h5adconttable -o h5ad_leiden_by_source.csv -x leiden -y source -- file.h5ad
	"""
	print('infile: ', infile)
	if not outfile:
		if bname:
			outfile = f"{bname}_{xvariable}_by_{yvariable}"
		else:
			outfile = f"h5ad_{xvariable}_by_{yvariable}"
	print('outfile: ', outfile)
	adata = sc.read(infile)
	df = pd.DataFrame(
		{
			xvariable:adata.obs[xvariable],
			yvariable:adata.obs[yvariable]
		}
	)
	table = pd.crosstab(index=df[xvariable], columns=df[yvariable])
	table.to_csv(f"{outfile}.csv", index=True)
	table.style.background_gradient(axis=None).to_excel(f"{outfile}.xlsx")

if __name__ == "__main__":
	main()
