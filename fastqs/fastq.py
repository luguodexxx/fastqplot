#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 4:05 PM
__author__ = 'Zhou Ran'

import sys
import logging
import gzip
from collections import defaultdict, Counter

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(asctime)s:%(message)s")


def gc(seq):
    """ Return the GC content of as an int
    >>> x = tuple('TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT')
    >>> gc(x)
    30
    """
    g = seq.count('G')
    c = seq.count('C')
    return int((g + c) / len(seq) * 100)


class FastqReader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """

    def __init__(self, f):

        if f.endswith('gz'):
            self.file = gzip.open(f)
            self._gz = True
        else:
            self.file = open(f)
            self._gz = False

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            # gzip file, if T then decode the string
            if self._gz:
                name = next(self.file).decode('utf-8').strip().split()[0]
                seq = next(self.file).decode('utf-8').strip()
                strand = next(self.file).decode('utf-8').strip()
                qual = next(self.file).decode('utf-8').strip()
            else:
                name = next(self.file).strip().split()[0]
                seq = next(self.file).strip()
                strand = next(self.file).strip()
                qual = next(self.file).strip()

            if name.count(':YM:Z:') > 0:
                tag, dtype, data = name.split(':')[-3:]
                name = ':'.join(name.split(':')[:-3])
                return Fastq(name=name, seq=seq, strand=strand, qual=qual, conv=data)
            else:
                return Fastq(name=name, seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n):
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                if name.count(':YM:Z:') > 0:
                    tag, dtype, data = name.split(':')[-3:]
                    name = ':'.join(name.split(':')[:-3])
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual, conv=data)
                else:
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual)

    def fileno(self):
        return self.file.fileno()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class Fastq(object):
    """
    A class to process the fastq line
    """

    def __init__(self, name='', seq='', strand='+', qual='', conv=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.conv = conv
        self.i = int()
        assert isinstance(name, str)
        assert isinstance(seq, str)
        assert isinstance(qual, str)

    def __iter__(self):
        return self

    def next(self):
        if self.i < len(self):
            value, self.i = self[self.i], self.i + 1
            return value
        else:
            raise StopIteration()

    def __getitem__(self, key):
        if self.conv:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key], self.conv[key])
        else:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key])

    def __next__(self):
        return self.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name[0] != '@':
            self.name = ''.join(['@', self.name])
        if self.conv:
            return '\n'.join(['{0}:YM:Z:{1}'.format(self.name, self.conv),
                              self.seq, self.strand, self.qual]) + '\n'
        else:
            return '\n'.join([self.name, self.seq, self.strand, self.qual]) + '\n'

    def __len__(self):
        return len(self.seq)

    def gc(self):
        """ Return the GC content of self as an int
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT', qual=' # # ##EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.gc()
        30
        """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return int((g + c) / len(self) * 100)


class Stats:
    """ Counter for characterization of NGS reads
    """

    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
        self.conv = defaultdict(lambda: defaultdict(int))

    def evaluate(self, seq, qual, conv=None):
        """ Evaluate read object at each position, and fill in nuc and qual dictionaries """
        self.gc[gc(seq)] += 1
        if conv:
            cpgs = cpg_map(seq)
        for i in range(1, len(seq) + 1):
            self.depth[i] += 1
            self.nuc[i][seq[i - 1]] += 1
            self.qual[i][qual[i - 1]] += 1
            if conv:
                if cpgs[i - 1] != 'N':
                    self.conv[i][conv[i - 1]] += 1

    def kmercount(self, seq, k=5):
        """ Count all kmers of k length in seq and update kmer counter.
        """
        for kmer in window(seq, n=k):
            self.kmers[kmer] += 1

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


if __name__ == '__main__':
    s = Stats()
    COUNTINFO = 1000000
    c_count = 1
    for line in FastqReader('R2.r1.fq'):
        if c_count % COUNTINFO == 0:
            logging.info('processed {} reads.'.format(c_count))
            c_count += 1

        s.evaluate(line.seq, line.qual)

    # print(list(map(lambda x: x["G"], s.nuc.values())))
    plot(s.nuc)
