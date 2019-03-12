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

REVERSE_DIC = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N"
}

NUC_ORDER = {
    "A": 1,
    "T": 2,
    "G": 3,
    "C": 4,
    "N": 5

}


def nucplot(countsdict,
            fileprefix=None,
            fig_bw={'figsize': (8, 6)},
            rev_axis=False,
            rev_nuc=False
            ):
    """

    :param countsdict:
    :param fileprefix:
    :param fig_bw:
    :param rev: rev means for reverse the axis numbers
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

        # here to reverse the base
        if rev_nuc:
            nuc = REVERSE_DIC[nuc]
        # here to add reverse the position
        if rev_axis:
            axes.plot(positions, [nuc_percent[pos][nuc] for pos in positions[::-1]])
        else:
            axes.plot(positions, [nuc_percent[pos][nuc] for pos in positions])

    axes.set_ylim(-5, 105)
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width, box.height])

    axes.set_axisbelow(True)
    axes.set_title('Base content (reversed)' if rev_nuc else 'Base content')
    axes.set_xlabel('Position (reversed)' if rev_axis else 'Position')
    axes.set_ylabel('Base content (% basecall)')

    # if rev_nuc:
    #     legendlabel = tuple((REVERSE_DIC[n] for n in nucs), key=lambda x: NUC_ORDER[x])
    #
    # else:
    #     legendlabel = tuple((n for n in nucs))

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
