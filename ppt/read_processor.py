#!/usr/bin/env python
# -*- coding: utf-8

class read_processor:
    def __init__(self, line):
        self.bc1 = 0
        self.bc2 = 0
        self.bc4 = 0
        self.bc = ""
        self.rawline = line
        self.sbc1 = ""
        self.sbc2 = ""
        self.sbc4 = ""
        self.bsbc4 = ""
        self.type = ""
        self.umi = ""
        self.dock = -1
        self.valid = False

    def align_score(self, str1, str2):
        score = 0
        for i in range(len(str1)):
            if str2[i] == 'N':
                continue
            if str2[i] == str1[i]:
                score += 2
            else:
                score -= 1
        return score

    def trim(self):
        if len(self.rawline) < 95:
            return
        self.dock = 0
        t = 0         # actually this is the shift for next barcode, not this one
        cur_s = 10
        bait = "GTGGCCGATGTTTCG"    # this is part of the linker between the first and second barcodes
        for i in range(-2, 7):
            qu = self.rawline[18 + i: 18 + i + 15]
            score = self.align_score(qu, bait)
            if score < cur_s:
                continue
            cur_s = score   # only keep the best score and its shift
            t = i
        if cur_s < 13:
            return
        self.dock = 1
        self.sbc1 = self.rawline[10: 18]  # this is the first barcode
        self.umi = self.rawline[:10]
        if t == -1:     # this means there is missing base at the beginning, so the first base of UMI is unknown
            self.umi = "N" + self.rawline[:9]
        elif t == -2:
            self.umi = "NN" + self.rawline[:8]

        # Second barcode
        bait = "ATCCACGTGCTTGAG"
        cur_s = 10
        tt = t
        for i in range(-2, 3):
            qu = self.rawline[56 + i + tt: 56 + i + tt + 15]
            score = self.align_score(qu, bait)
            if score < cur_s:
                continue
            cur_s = score
            tt = t + i
        if cur_s < 12:
            return
        t = tt
        self.sbc2 = self.rawline[48 + t: 56 + t]
        self.dock = 2

        # Fourth barcode
        bait = "AGGCCAGAGCATTCG"
        cur_s = 10
        for i in range(-2, 3):
            qu = self.rawline[71 + i + tt: 71 + i + tt + 15]
            score = self.align_score(qu, bait)
            if score < cur_s:
                continue
            cur_s = score
            tt = t + i
        if cur_s < 11:
            return
        t = tt
        if len(self.rawline) < 89 + t:
            return
        self.sbc4 = self.rawline[87 + t: 91 + t]
        self.dock = 4

        # Type determination
        #  DNA AG
		#  RNA TC
        cur_s = 0
        b1 = self.rawline[86 + t]
        b2 = self.rawline[91 + t]
        if b1 == "A":
            cur_s += 1
        if b1 == "T":
            cur_s -= 1
        if b2 == "G":
            cur_s += 1
        if b2 == "C":
            cur_s -= 1
        if cur_s > 0:         
            self.type = "d"
        elif cur_s < 0:
            self.type = "r"
        else:
            self.type = "n"    

    