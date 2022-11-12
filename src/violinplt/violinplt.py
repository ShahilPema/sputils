#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import click
import os
from pathlib import Path
from exeR.exeR import exeR

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('infile', type=click.Path(exists=True), nargs=1)
@click.option('-o', '--outdir', type=click.Path(exists=True), required=True, help='The output directory for violin plots')
@click.option('-b', '--bname', type=click.Path(), required=True, help='Bname.')
@click.option('-v', '--verbose', is_flag=True, help='Verbose.')
def main(infile, outdir, bname, verbose):
	"""
Create a violin plot for given samples.

\b
Example:
	violinplt.py sample.h5seurat -o .

\b
Date: 2022/11/04
Authors: Shahil Pema <pemashahil@yahoo.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	filename=Path(__file__).stem
	script=f'{absdir}/R/{filename}.R'
	exprs=[
	f"infile='{infile}'",
	f"outdir='{outdir}'",
	f"bname='{bname}_dfviolin.pdf'"
	]
	return exeR.callback(exprs, script=script, condaenv='', verbose=True)

if __name__ == "__main__":
	sys.exit(main())
