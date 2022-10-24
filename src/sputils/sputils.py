#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import click
from importlib_metadata import packages_distributions

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-l', '--listets', is_flag=True, help='List entry points (i.e., utilities). (default: False)')
@click.version_option()
def main(listets):
	"""
Wrappers in Python.

\b
Example:
  sputils -l
	"""
	if listets:
		for k,v in packages_distributions().items():
			if 'sputils' in v:
				print(','.join(v), k, sep='\t')

if __name__ == "__main__":
	main()
