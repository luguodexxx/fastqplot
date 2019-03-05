#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 5:13 PM
__author__ = 'Zhou Ran'

import sys
import click

from .fastq import Stats, FastqReader
from .plot import nucplot


@click.command()
@click.option('--fq',
              type=str,
              help='The fastq file input, gzip or not alse fine')
@click.option('--fp',
              type=str,
              help='The outputfile prefix(fp).')
@click.option('--rev',
              is_flag=True,
              help='Reverse the base location, defaule:False.'
              )
def plot(fq, fp, rev):
    """Stat the fastq file and plot the base distribution plot"""

    if not all([fq, fp]):
        plot(['plot', '--help'])
        sys.exit(1)

    s = Stats()
    for line in FastqReader(fq):
        s.evaluate(line.seq, line.qual)

    nucplot(s.nuc,
            fileprefix=fp,
            rev=rev)


if __name__ == '__main__':
    plot()
