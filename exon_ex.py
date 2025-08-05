#!/usr/bin/env python3

from sys import stderr, exit
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType


def parse_attributes_gff(attributes):
    """Parse GFF3 attributes like ID=...;Parent=..."""
    attr_dict = {}
    for field in attributes.strip().split(';'):
        if '=' not in field:
            continue
        key, value = field.strip().split('=', 1)
        attr_dict[key] = value
    return attr_dict


def extract_exons(gtf_file, verbose=False):
    genes = dd(list)
    trans = {}

    # Parse valid exon lines from the GFF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature.lower() != 'exon' or left >= right:
            continue

        values_dict = parse_attributes_gff(values)

        # Accept either transcript_id or Parent (common in GFF3)
        transcript_id = values_dict.get('transcript_id') or values_dict.get('Parent')
        gene_id = values_dict.get('gene_id') or values_dict.get('gene') or transcript_id

        if not transcript_id or not gene_id:
            continue

        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[gene_id].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    tmp_exons = set()
    for chrom, strand, texons in trans.values():
        for i in range(len(texons)):
            tmp_exons.add((chrom, texons[i][0], texons[i][1], strand))
    tmp_exons = sorted(tmp_exons)
    if len(tmp_exons) <= 0:
        return

    exons = [tmp_exons[0]]
    for exon in tmp_exons[1:]:
        prev_exon = exons[-1]
        if exon[0] != prev_exon[0]:
            exons.append(exon)
            continue
        assert prev_exon[1] <= exon[1]
        if prev_exon[2] < exon[1]:
            exons.append(exon)
            continue

        if prev_exon[2] < exon[2]:
            strand = prev_exon[3]
            if strand not in "+-":
                strand = exon[3]
            exons[-1] = (prev_exon[0], prev_exon[1], exon[2], strand)

    for chrom, left, right, strand in exons:
        # Zero-based offset
        print('{}\t{}\t{}\t{}'.format(chrom, left - 1, right - 1, strand))


if __name__ == '__main__':
    parser = ArgumentParser(description='Extract exons from a GFF file')
    parser.add_argument('gtf_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GFF file (use "-" for stdin)')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    extract_exons(args.gtf_file, args.verbose)

