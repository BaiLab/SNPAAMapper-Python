#!/usr/bin/env python
'''
Algorithm for mapping identified variants onto the genomic location and
reporting the hit class

USAGE:
    python algorithm_mapping_variants_reporting_class_intronlocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf
'''

import sys, csv
from sys import argv

EXON_FILE = None
SNP_FILE = None
EXON_BUFFER = None
INTRON_OPTION = None

# increase field limit
LIMIT = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(LIMIT)
        break
    except OverflowError: LIMIT = int(LIMIT/10)


if len(argv) == 3:
    EXON_FILE = argv[1]
    SNP_FILE = argv[2]
    print('The program assumes that you do NOT want to report how far the variant falls in the exon boundary.')
elif len(argv) == 4:
    EXON_FILE = argv[1]
    SNP_FILE = argv[2]
    EXON_BUFFER = int(argv[3])
    print('The program assumes that you DO want to report how far the variant falls in the exon boundary.')
    print('Only variants flanking their nearby exon within <=' + str(EXON_BUFFER) + ' bp are reported')
    INTRON_OPTION = True
else:
    print('The input commands do not meet the requirement. Please see the README file and try again.')
    sys.exit()

output_file = SNP_FILE + '.append'

class FeatureType:
    """class for reading in exon files"""
    def __init__(self, name, ext):
        self.name = name
        self.ext = ext
        self.chrom = {} # start position -> stop position
        self.gene = {} # chrom_start -> stop position
        self.starts = {} # chrom -> start positions

        with open(EXON_FILE + '.' + ext + '_link_shrink', encoding='utf-8') as linkshr:
            for lin in list(csv.reader(linkshr, delimiter='\t')):
                self.chrom[lin[0]] = int(lin[1])
        with open(EXON_FILE + '.' + ext + '_gene', encoding='utf-8') as gen:
            for lin in list(csv.reader(gen, delimiter='\t')):
                self.gene[lin[0]] = lin[1]

        # store start positions for binary search
        with open(EXON_FILE + '.' + ext, encoding='utf-8') as exon:
            for lin in list(csv.reader(exon, delimiter='\t')):
                chrom = lin[0]
                arr = []
                for start in lin[1:]:
                    if start != '': arr.append(int(start))
                if chrom not in self.starts: self.starts[chrom] = arr
                else: self.starts.extend(arr)


def binary_search(start, snp_st):
    '''faster binary search'''
    l_o = 0
    h_i = len(start) - 1
    while l_o < h_i:
        midd = (l_o + h_i + 1) // 2
        if snp_st >= start[midd]: l_o = midd
        else: h_i = midd - 1
    return l_o

feature_types = [FeatureType('CDSHIT', 'cds'), FeatureType('INTRON', 'intron'),
                 FeatureType('UTR5', 'utr5'), FeatureType('UTR3', 'utr3'),
                 FeatureType('UPSTREAM', 'upstream'),
                 FeatureType('DOWNSTREAM', 'downstream')]

# CDS - start is 0-indexed, end is 1-indexed
# UTR5 - start is 0-indexed, end is 0-indexed
# UPSTREAM - start is 0-indexed, end is 0-indexed
# UTR3 - start is 1-indexed, end is 1-indexed
# DOWNSTREAM - start is 1-indexed, end is 1-indexed
# INTRON - start is 1-indexed, add buffer because end is 0-indexed

with open(output_file, 'w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file, delimiter='\t')

    # map vcf variant back onto genome location and report hit class
    with open(SNP_FILE, encoding='utf-8') as snp:
        for line in list(csv.reader(snp, delimiter='\t')):
            if line[0][0] == '#': continue
            SNP_CHR = line[0]
            if 'chr' not in SNP_CHR: SNP_CHR = 'chr' + SNP_CHR
            if SNP_CHR == 'chrMT': SNP_CHR = 'chrM'
            snp_start = int(line[1])

            for typ in feature_types:
                starts = typ.starts[SNP_CHR]
                mid = binary_search(starts, snp_start-1)
                output = [SNP_CHR,*line[1:].copy()]

                if typ.name == 'INTRON' and INTRON_OPTION:
                    if starts[mid] <= snp_start < starts[mid] + EXON_BUFFER or typ.chrom[SNP_CHR + '_' + str(starts[mid])] + 1 - EXON_BUFFER < snp_start <= typ.chrom[SNP_CHR + '_' + str(starts[mid])] + 1:
                        if snp_start - starts[mid] < EXON_BUFFER:
                            output.append('INTRONHIT.' + str(snp_start - starts[mid] + 1) + '.right')
                            output.append(typ.gene[SNP_CHR + '_' + str(starts[mid])])
                            writer.writerow(output)
                        elif feature_types[1].chrom[SNP_CHR + '_' + str(starts[mid])] + 1 - snp_start < EXON_BUFFER:
                            output.append('INTRONHIT.' + str(typ.chrom[SNP_CHR + '_' + str(starts[mid])] - snp_start + 2))
                            output.append(typ.gene[SNP_CHR + '_' + str(starts[mid])])
                            writer.writerow(output)
                elif snp_start-1 >= starts[mid] and snp_start <= typ.chrom[SNP_CHR + '_' + str(starts[mid])]:
                    output.append(typ.name)
                    output.append(typ.gene[SNP_CHR + '_' + str(starts[mid])])
                    writer.writerow(output)

            print("Done for SNP " + SNP_CHR + "_" + str(snp_start) + ".")
