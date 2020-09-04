#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 5:13 PM
__author__ = 'Zhou Ran'

import sys
import click
import logging
import pandas as pd

from .fastq import Stats, FastqReader, statnucfromfile
from .plot import nucplot

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(asctime)s:%(message)s")


@click.group()
def cli():
    """A command line tools to plot fastq stats"""
    pass


@click.command()
@click.option('--fq',
              type=str,
              help='The fastq file input, gzip or not alse fine')
@click.option('--fp',
              type=str,
              help='The outputfile prefix(fp).')
@click.option('--maxlen',
              type=int,
              default=300,
              help='The max length of reads for stat')

@click.option('--revx',
              is_flag=True,
              help='Reverse the base location(x axis), defaule:False.'
              )
@click.option('--revn',
              is_flag=True,
              help='Reverse the nucleobase, defaule:False.'
              )
def fqplot(fq,
           fp,
           maxlen,
           revx,
           revn):
    """Stat the fastq file and plot the base distribution plot"""

    if not all([fq, fp]):
        cli(['fqplot', '--help'])
        sys.exit(1)

    COUNTINFO = 1000000
    c_count = 1

    s = Stats()
    for line in FastqReader(fq):
        if c_count % COUNTINFO == 0:
            logging.info('processed {} reads.'.format(c_count))
        c_count += 1
        s.evaluate(line.seq, line.qual, max_len=maxlen)

    pdm = pd.DataFrame.from_dict(s.nuc)
    pdm.to_csv('{}.nulc_dis.txt'.format(fp), na_rep=0)
    nucplot(s.nuc,
            fileprefix=fp,
            rev_axis=revx,
            rev_nuc=revn)


@click.command()
@click.option('--file',
              type=str,
              help='The nuc config file')
@click.option('--fp',
              type=str,
              help='The outputfile prefix(fp).')
@click.option('--revx',
              is_flag=True,
              help='Reverse the base location(x axis), defaule:False.'
              )
@click.option('--revn',
              is_flag=True,
              help='Reverse the nucleobase, defaule:False.'
              )
def nucfile(
        file,
        fp,
        revx,
        revn):
    """Fast way to plot the saved information"""
    if not all([file, fp]):
        cli(['nucfile', '--help'])
        sys.exit(1)

    s = statnucfromfile(file)
    nucplot(s,
            fileprefix=fp,
            rev_axis=revx,
            rev_nuc=revn)


cli.add_command(fqplot)
cli.add_command(nucfile)

if __name__ == '__main__':
    cli()
