#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 5:10 PM
__author__ = 'Zhou Ran'

import sys
from collections import defaultdict
import matplotlib as mpl

if sys.platform is not 'darwin':
    mpl.use('Agg')
import matplotlib.pyplot as plt


def nucplot(countsdict,
            fileprefix=None,
            fig_bw={'figsize': (8, 6)},
            rev=False):
    """

    :param countdict:
    :return:
    """
    nucs = [
        "A",
        "T",
        "G",
        "C",
        "N"
    ]
    legend_bg_color_kw = 'facecolor'
    fig, axes = plt.subplots(
        nrows=1, subplot_kw={legend_bg_color_kw: 'white'}, **fig_bw)

    nuc_percent = defaultdict(lambda: defaultdict(int))

    positions = list(countsdict.keys())

    for pos, count in tuple(countsdict.items()):
        max_depth = sum(tuple(count.values()))
        for nuc in nucs:
            if max_depth > 0:
                nuc_percent[pos][nuc] = float(count[nuc]) / max_depth * 100
            else:
                nuc_percent[pos][nuc] = 0.
    for nuc in nucs:
        if rev:
            axes.plot(positions, [nuc_percent[pos][nuc] for pos in positions[::-1]])
        else:
            axes.plot(positions, [nuc_percent[pos][nuc] for pos in positions])

    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width, box.height])

    axes.set_axisbelow(True)
    axes.set_title('Base content')
    axes.set_xlabel('Position')
    axes.set_ylabel('Base content (% basecall)')

    legend = axes.legend(
        tuple((n for n in nucs)),
        ncol=len(nucs),
        bbox_to_anchor=(1, 1),
        loc=1,  # for best
        prop={
            'size': 8
        },
        frameon=False
    )
    frame = legend.get_frame()
    frame.set_facecolor('white')
    for label in legend.get_texts():
        label.set_color('black')

    if fileprefix:
        plt.savefig('{}.basedistribution.pdf'.format(fileprefix))
    else:
        plt.savefig('BaseDistribution.pdf')
