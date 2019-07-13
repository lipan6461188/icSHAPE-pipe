#-*- coding:utf-8 -*-

import re, pysam


def reverse_comp(seq):
    reverseDict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N',
               'a':'t', 'c':'g', 't':'a', 'g':'c', '-':'-'}
    return ''.join(map(lambda x: reverseDict[x], list(seq[::-1])))

def cutSeq(rawSeq, lineLen=60):
    cut_seq = ''
    idx = 0
    while idx < len(rawSeq):
        cut_seq += rawSeq[idx:idx+lineLen]+'\n'
        idx += lineLen
    return cut_seq[:-1]

class SeqClass(object):
    def __init__(self, fa_fileName):
        self.fileName = fa_fileName
        self.genome = pysam.Fastafile(fa_fileName)
        print("seqClass: 输入坐标前闭后开，0-based")
    def fetch(self, Chr, Start, End, Strand="+"):
        """获取序列"""
        if Strand == '+':
            return self.genome.fetch(Chr, Start, End)
        if Strand == '-':
            return reverse_comp(self.genome.fetch(Chr, Start, End))
    def has(self, Chr):
        return Chr in self.genome.references



