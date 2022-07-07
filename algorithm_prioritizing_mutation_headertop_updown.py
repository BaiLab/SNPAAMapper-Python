#!/usr/bin/env python
'''
Algorithm for prioritizing mutation effects

USAGE:
    python algorithm_prioritizing_mutation_headertop_updown.py VCF_input_file_in_tab_delimited_format.txt.append.out.txt
'''

import sys, csv
from sys import argv

input_file = argv[1]
output_file = input_file + ".prioritized_out"

# increase field limit
LIMIT = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(LIMIT)
        break
    except OverflowError: LIMIT = int(LIMIT/10)


with open(input_file, mode='r', encoding='utf-8') as line: line_array = line

with open(input_file, encoding='utf-8') as vcf:
    data = list(csv.reader(vcf, delimiter='\t'))
columns = data[0]
data = data[1:]

score = []
output_row = {} # stores highest-scoring row number for each variant

for i, line in enumerate(data):
    score.append(-1)
    if line[12] == 'CDSHIT':
        # get rid of right reciprocal pairs
        if line[7] == 'SNP':
            for s in line[9].split(','):
                if s == 'NSM': score[i] = 6
                if s == 'NSN': score[i] = 5
                if s == 'SYN': score[i] = 4
        else: score[i] = 3.5
    elif line[12] == 'UPSTREAMHIT': score[i] = 3
    elif line[12] == 'UTR3HIT': score[i] = 2
    elif line[12] == 'UTR5HIT': score[i] = 1
    elif line[12] == 'DOWNSTREAMHIT': score[i] = 0

    # use combination key: sample_chromosome_gene
    KEY = "_".join(line[0:3])

    # save row number with the highest score
    if KEY not in output_row or score[i] > score[output_row[KEY]]:
        output_row[KEY] = i

with open(output_file, 'w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(columns)
    for i, row in output_row.items(): writer.writerow(data[row])
