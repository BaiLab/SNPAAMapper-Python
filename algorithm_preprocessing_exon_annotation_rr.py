#!/usr/bin/env python
'''
Algorithm for preprocessing the gene structure to build annotation structure for
each exon

USAGE:
    python algorithm_preprocessing_exon_annotation_rr.py ChrAll_knownGene.txt.exon
'''

from sys import argv
import csv

class FeatureType:
    """ Feature Type """

    def __init__(self, ext):
        self.ext = ext
        self.file = open(exon_file + '.' + ext, 'w', newline='', encoding='utf-8')
        self.gene_file = open(exon_file + '.' + ext + '_gene', 'w', newline='', encoding='utf-8')
        self.link_file = open(exon_file + '.' + ext + '_link', 'w', newline='', encoding='utf-8')
        self.link_shrink_file = open(exon_file + '.' + ext + '_link_shrink', 'w', newline='', encoding='utf-8')
        self.chrom_starts = {}
        self.chrom_ends = {}
        self.chrom_genes = {}
        self.sorted_chrom = []

    def process_exon(self, start, end, chrm, exon_gen):
        """ append start position, end position, and gene to lists that chrom is
        mapped to """
        if start != 'NA' or end != 'NA':
            if chrm in self.chrom_starts:
                self.chrom_starts[chrm].append(int(start))
                self.chrom_ends[chrm].append(int(end))
                self.chrom_genes[chrm].append(exon_gen)
            else:
                self.chrom_starts[chrm] = [int(start)]
                self.chrom_ends[chrm] = [int(end)]
                self.chrom_genes[chrm] = [exon_gen]

    def output_info(self):
        """ print start/end position and gene to _link file """
        for chrm, cont in self.chrom_starts.items():
            writer = csv.writer(self.link_file, delimiter='\t')
            for i, lin in enumerate(cont):
                writer.writerow([chrm + '_' + str(lin),
                                 self.chrom_ends[chrm][i],
                                 self.chrom_genes[chrm][i]])
        self.link_file.close()

    def output_sorted(self):
        """ output chromosomes and start positions in sorted order """
        writer = csv.writer(self.file, delimiter='\t')
        for key in sorted(self.chrom_starts):
            sort_chrom_array = sorted(self.chrom_starts[key])
            output = [key]
            output.extend(sort_chrom_array)
            writer.writerow(output)
        self.file.close()

    def output_max(self):
        """ output """
        with open(exon_file + '.' + self.ext + '_link', encoding='utf-8') as ex:
            exons = list(csv.reader(ex, delimiter='\t'))
        output_row = {} # start position -> row with greatest stop position
        link_shrink_writer = csv.writer(self.link_shrink_file, delimiter='\t')
        gene_writer = csv.writer(self.gene_file, delimiter='\t')

        for i, lin in enumerate(exons):
            start = lin[0]
            stop = lin[1]
            #gene = lin[2]

            if start not in output_row or int(stop) > int(exons[output_row[start]][1]):
                output_row[start] = i

        for start, row in output_row.items():
            link_shrink_writer.writerow([start, exons[row][1]])
            gene_writer.writerow([start, exons[row][2]])

        self.link_shrink_file.close()
        self.gene_file.close()


exon_file = argv[1]
feature_types = [FeatureType('cds'), FeatureType('intron'), FeatureType('utr5'),
                 FeatureType('utr3'), FeatureType('upstream'),
                 FeatureType('downstream')]

with open(exon_file, encoding='utf-8') as exo:
    for exon_line in list(csv.reader(exo, delimiter = '\t')):
        chrom = exon_line[2]
        if exon_line[0] == 'bin' or '_' in chrom or 'chrom' in chrom: # ignore header and strange chromosome pieces
            continue
        exon_gene = exon_line[1]
        feature_types[0].process_exon(exon_line[5], exon_line[6], chrom, exon_gene) # cds
        feature_types[1].process_exon(exon_line[10], exon_line[11], chrom, exon_gene) # intron
        feature_types[2].process_exon(exon_line[12], exon_line[13], chrom, exon_gene) # utr5
        feature_types[3].process_exon(exon_line[14], exon_line[15], chrom, exon_gene) # utr3
        feature_types[4].process_exon(exon_line[16], exon_line[17], chrom, exon_gene) # upstream
        feature_types[5].process_exon(exon_line[18], exon_line[19], chrom, exon_gene) # downstream

for typ in feature_types:
    typ.output_info()
    typ.output_sorted()
    typ.output_max()
