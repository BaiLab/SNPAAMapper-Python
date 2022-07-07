#!/usr/bin/env python
'''
Algorithm for predicting amino acid changes

USAGE:
    python algorithm_predicting_full_aa_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.vcf.out.txt
'''

import sys, csv, re
from sys import argv

snp_file = argv[1]
conversion_file = argv[2]
gene_outfile = argv[3]
cds_intron_file = argv[4]
output_file = snp_file + '.out.txt'

# increase field limit
LIMIT = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(LIMIT)
        break
    except OverflowError: LIMIT = int(LIMIT/10)

# codon -> AA conversion
aa_dict = {'ATG':'M', 'TGG':'W', 'TTT':'F', 'TTC':'F', 'TAT':'Y', 'TAC':'Y',
           'TGT':'C', 'TGC':'C', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
           'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D',
           'GAA':'E', 'GAG':'E', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'CCT':'P',
           'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T',
           'ACG':'T', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A',
           'GCC':'A', 'GCA':'A', 'GCG':'A', 'GGT':'G', 'GGC':'G', 'GGA':'G',
           'GGG':'G', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S',
           'AGC':'S', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L',
           'CTG':'L', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R',
           'AGG':'R', 'TAG':'*', 'TAA':'*', 'TGA':'*',}

# UCSC_ID -> gene symbol
gene_dict = {}
# UCSC_ID -> strand
strand_dict = {}
# UCSC_ID -> start, end, strand, cds line
cds_dict = {}

with open(conversion_file, encoding='utf-8') as conversions:
    for line in list(csv.reader(conversions, delimiter = '\t')):
        gene_dict[line[0]] = line[4]
with open(cds_intron_file, encoding='utf-8') as cds_introns:
    for line in list(csv.reader(cds_introns, delimiter = '\t')):
        strand_dict[line[0]] = line[2]

with open(gene_outfile, encoding='utf-8') as genes:
    # HANDLE ONLY single-line-fashion sequences in hg19_CDSIntronWithSign.txt.out ------------------#
    # IF single-line-fashion, the following line is faster than the method that handles both fashions
    #gene_out_list = list(csv.reader(genes, delimiter = ' '))
    # HANDLE BOTH single- and multiple-line-fashion sequences in hg19_CDSIntronWithSign.txt.out ----#
    sequences = genes.read()
    sequences = re.split("^>", sequences, flags=re.MULTILINE)
    del sequences[0]
gene_out_list=[]
for fasta in sequences:
    try: header, sequence = fasta.split("\n", 1)
    except ValueError: print(fasta)
    header = ">" + header
    sequence = sequence.replace("\n","")
    gene_out_list.append(header.split())
    if sequence: gene_out_list.append(sequence.split())
    #-----------------------------------------------------------------------------------------------#

for i in range(0, len(gene_out_list), 2):
    if len(gene_out_list[i]) == 6:
        ucsc_id = gene_out_list[i][0][16:]
        strand = gene_out_list[i][4][-1]
        arr_chr = gene_out_list[i][1].split(':')[0].split('=')
        arr = gene_out_list[i][1].split(':')[1].split('-')
        cds_start = int(arr[0])
        cds_end = int(arr[1])

        if ucsc_id not in cds_dict: cds_dict[ucsc_id] = {}
        cds_dict[ucsc_id][arr_chr[1]] = (cds_start, cds_end, strand, gene_out_list[i+1][0])


def protein_translation(original_cds_line, cds_lin, snp_lo, protein_flg, strad, case_ct, cases_left, out):
    '''function for translating protein'''
    codon_change_str = ''
    truncate_original_cds_line = ''
    truncate_original_aa_line = ''
    checked_codon = None

    # translate CDS into AA sequence by moving 3 at a time
    for j in range(0, len(original_cds_line), 3):
        original_codon = original_cds_line[j:j+3]
        truncate_original_cds_line = truncate_original_cds_line + original_codon
        truncate_original_aa_line = truncate_original_aa_line + aa_dict[original_codon]

        # if this is the codon position
        if snp_lo - j >= 0 and snp_lo - j <= 2 and protein_flg:
            checked_codon = aa_dict[original_codon]
            # first case of ALT
            if case_ct in (1, cases_left):
                out.extend([strad, j//3 + 1, 'SNP', checked_codon + '(' + original_codon + ')'])
            else: out.extend([checked_codon + '(' + original_codon + ')'])

    truncate_cds_line = ''
    truncate_aa_line = ''
    for j in range(0, len(cds_lin), 3):
        codon = cds_lin[j:j+3]
        truncate_cds_line = truncate_cds_line + codon
        try:
            truncate_aa_line = truncate_aa_line + aa_dict[codon]

            if snp_lo - j >= 0 and snp_lo - j <= 2 and protein_flg:
                #snp_codon_mutant = codon
                snp_aa_mutant = aa_dict[codon]

                out[-1] = out[-1] + '->' + snp_aa_mutant + '(' + codon + ')'
                if case_ct == 1:
                    if checked_codon == snp_aa_mutant: out.append('SYN')
                    elif snp_aa_mutant == '*': out.append('NSN')
                    else: out.append('NSM')
                else:
                    out[-1] = out[-1] + ','
                    if checked_codon == snp_aa_mutant:
                        codon_change_str = codon_change_str + 'SYN'
                    elif snp_aa_mutant == '*':
                        codon_change_str = codon_change_str + 'NSN'
                    else: codon_change_str = codon_change_str + 'NSM'
                    if cases_left > 1:  # more ALT cases left
                        codon_change_str = codon_change_str  + ','
        except KeyError: print(codon+"\tinvldCodon")

    print(truncate_original_cds_line)
    print(truncate_cds_line)
    print(truncate_original_aa_line)
    print(truncate_aa_line)

    # newly added for printing full length of AA
    out.append(truncate_original_aa_line)
    out.append(truncate_aa_line)

    return codon_change_str

def remove_intron(lin):
    '''function for removing intron in line'''
    return re.sub(r'[acgtn]', '', lin)

def rev_dna_comp(dna):
    '''function for translating'''
    return dna[::-1].translate(dna.maketrans('ACGTacgt', 'TGCAtgca'))

with open(output_file, 'w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Sample', 'Chromosome', 'Variant Position', 'Gene Symbol',
                     'UCSC ID', 'Strand', 'AA Position of Mutation (for CDSHIT)',
                     'Variant Type', 'Amino Acid Ref (Codon) -> AA SNP (Codon)',
                     'Variant Class', 'Ref AA chain', 'Alt AA chain', 'Hit Type',
                     'Known dbSNP', 'Ref nt', 'Alt nt', 'Quality', 'Depth',
                     'Allele Freq', 'Read Categories', 'Info'])

    # read SNP file
    with open(snp_file, encoding='utf-8') as snps:
        reader=list(csv.reader(snps, delimiter = '\t'))
    # in case there is a header line in the snp_file (e.g. 007_crop.vcf.append)
    if "CHROM" in list(reader)[0][0]: reader1=reader[1:]
    else: reader1=reader
    for line in reader1:
        SNP_FLAG = True  # initially assume line is an SNP
        PROTEIN_FLAG = False

        snp_chromosome = line[0]
        snp_loc = int(line[1])

        # assume ref has only one case
        ref_char = line[3].upper()
        snp_chars = line[4].upper().split(',')
        case_count = len(snp_chars)
        for s in snp_chars:
            if len(s) > 1:  # not a single SNP
                SNP_FLAG = False
                break

        # get depth and read category info
        info_array = line[7].split(';')
        depth_array = []
        alle_freq = []
        read_category = []

        HIT_TYPE = None
        UCSC_ID = None

        if 'VDB' in line[7]:  # samtools v0.1.18
            if len(line) == 14:  # 3 samples run from samtools
                HIT_TYPE = line[12]
                UCSC_ID = line[13]
                print('3 samples, old samtools..')
            elif len(line) == 13:  # 2 samples run from samtools
                HIT_TYPE = line[11]
                UCSC_ID = line[12]
                print('2 samples, old samtools..')
            elif len(line) == 12:  # 1 sample run from samtools
                HIT_TYPE = line[10]
                UCSC_ID = line[11]
                print('1 sample, old samtools..')
        else: # e.g. clinvar data
            print('1 sample, new samtools..')
            HIT_TYPE = line[8]
            UCSC_ID = line[9]

        UCSC_ID_FLAG = False

        STRAND = None
        CODON_CHANGE_STRING = ''
        output = []
        if HIT_TYPE == 'CDSHIT':
            print('\n===========================')
            print(*line, gene_dict[UCSC_ID], sep='\t')
            print('The number of possible ALT cases for this line is ' + str(case_count))
            # both ref and SNP calls are single SNPs
            if len(ref_char) == 1 and SNP_FLAG: PROTEIN_FLAG = True
            output = [snp_file.split('.')[0], line[0][3:], line[1], gene_dict[UCSC_ID], UCSC_ID]

            # for each possible ALT case
            for case in range(0, case_count):
                if case > 0: print('\n--------------------------------------')
                CDS_START = None
                CDS_END = None
                LINE_FLAG = False

                if UCSC_ID in cds_dict:
                    if snp_chromosome in cds_dict[UCSC_ID]:
                        if snp_loc < cds_dict[UCSC_ID][snp_chromosome][0] or snp_loc > cds_dict[UCSC_ID][snp_chromosome][1]:
                            print('SNP location is outside of CDS region! No checking occurs!')

                        CDS_START = cds_dict[UCSC_ID][snp_chromosome][0]
                        CDS_END = cds_dict[UCSC_ID][snp_chromosome][1]
                        STRAND = cds_dict[UCSC_ID][snp_chromosome][2]
                        line2 = cds_dict[UCSC_ID][snp_chromosome][3]
                        LINE_FLAG = True
                        UCSC_ID_FLAG = True

                # found CDS line
                if LINE_FLAG:
                    cds_line = line2
                    line2_no_intron = remove_intron(line2)
                    BEFORE_SNP_STRING = None

                    # replace original char with SNP
                    if STRAND == '+':
                        coordinate = snp_loc - CDS_START
                        # count number of lowercase nucleotides (NON-CDS or intron nucleotide) before SNP location and adjust its corrdinate
                        BEFORE_SNP_STRING = line2[:coordinate]
                        cds_line = cds_line[:coordinate] + snp_chars[case] + cds_line[coordinate + len(ref_char):]

                        # process new CDS string
                        cds_line_no_intron = remove_intron(cds_line)
                        if case_count == 1:
                            protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(re.findall(r'[a-z]', BEFORE_SNP_STRING)), PROTEIN_FLAG, STRAND, case_count, case_count - case, output)
                        else:
                            CODON_CHANGE_STRING = CODON_CHANGE_STRING + protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(re.findall(r'[a-z]', BEFORE_SNP_STRING)), PROTEIN_FLAG, STRAND, case_count, case_count - case, output)
                    else:
                        coordinate = CDS_END - snp_loc
                        replaced_snp_char = rev_dna_comp(snp_chars[case]).upper()
                        BEFORE_SNP_STRING = line2[:coordinate]

                        # process new CDS string
                        cds_line = cds_line[:coordinate - len(ref_char) + 1] + replaced_snp_char + cds_line[coordinate + 1:]
                        cds_line_no_intron = remove_intron(cds_line)
                        if case_count == 1:
                            protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(ref_char) + 1 - len(re.findall(r'[a-z]', BEFORE_SNP_STRING)), PROTEIN_FLAG, STRAND, case_count, case_count - case, output)
                        else:
                            CODON_CHANGE_STRING = CODON_CHANGE_STRING + protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(ref_char) + 1 - len(re.findall(r'[a-z]', BEFORE_SNP_STRING)), PROTEIN_FLAG, STRAND, case_count, case_count - case, output)
                    LINE_FLAG = False
                    break

            if len(ref_char) == 1 and SNP_FLAG:  # single SNP case
                depth_array = info_array[0].split('=')
                if 'VDB' in line[7]: # samtools v0.1.18
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                else: # samtools v0.1.12
                    alle_freq = info_array[1].split('=')
                    read_category = info_array[3].split('=')
                if len(CODON_CHANGE_STRING) == 0:  # no ALT cases
                    output.extend([HIT_TYPE, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
                else:
                    output.extend([CODON_CHANGE_STRING, HIT_TYPE, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            else: # indel case
                depth_array = info_array[1].split('=')
                if 'VDB' in line[7]:  # samtools v0.1.18
                    alle_freq = info_array[3].split('=')
                    read_category = info_array[5].split('=')
                else:  # samtools v0.1.12
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                output.extend([STRAND, '---', 'INDEL', '---', '---', '---', '---', HIT_TYPE, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            if not UCSC_ID_FLAG: print(UCSC_ID + ' cannot be found!!')
        else:  # INTRON, UTR, or UP-DOWNSTREAM
            output.extend([snp_file.split('.')[0], line[0][3:], line[1], gene_dict.get(UCSC_ID, 'Missing'), UCSC_ID])

            if len(ref_char) == 1 and SNP_FLAG:  # both ref and SNP calls are single SNPs
                depth_array = info_array[0].split('=')
                if 'VDB' in line[7]:  # samtools v0.1.18
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                else:  # samtools v0.1.12
                    alle_freq = info_array[1].split('=')
                    read_category = info_array[3].split('=')
                output.extend([strand_dict[UCSC_ID], '---', 'SNP', '---', '---', '---', '---', HIT_TYPE, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            else:
                depth_array = info_array[1].split('=')
                if 'VDB' in line[7]:  # samtools v0.1.18
                    alle_freq = info_array[3].split('=')
                    read_category = info_array[5].split('=')
                else:  # samtools v0.1.12
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                output.extend([strand_dict[UCSC_ID], '---', 'INDEL', '---', '---', '---', '---', HIT_TYPE, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
        writer.writerow(output)
