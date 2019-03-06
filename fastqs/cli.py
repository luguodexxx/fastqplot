#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 5:13 PM
__author__ = 'Zhou Ran'

import sys
import click

from .fastq import Stats, FastqReader
from .plot import nucplot

if __name__ == '__main__':
    plot()


@click.command()
@click.option('--fq',
              type=str,
              help='The fastq file input, gzip or not alse fine')
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
def plot(fq,
         fp,
         revx,
         revn):
    """Stat the fastq file and plot the base distribution plot"""

    if not all([fq, fp]):
        plot(['plot', '--help'])
        sys.exit(1)

    COUNTINFO = 1000000
    c_count = 1

    s = Stats()
    for line in FastqReader(fq):
        if c_count % COUNTINFO == 0:
            logging.info('processed {} reads.'.format(c_count))
            c_count += 1
        s.evaluate(line.seq, line.qual)

    # print(list(map(lambda x: x["G"], s.nuc.values())))

    nucplot(s.nuc,
            fileprefix=fp,
            rev_axis=revx,
            rev_nuc=revn)
