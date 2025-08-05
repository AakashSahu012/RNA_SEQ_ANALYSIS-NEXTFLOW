#!/usr/bin/env python3

from sys import stderr, exit
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType

def parse_attributes_gff3(attr_str):
    attr_dict = {}
    for field in attr_str.strip().split(';'):
        if '=' in field:
            key, value = field.strip().split('=', 1)
            attr_dict[key] = value
    return attr_dict

def extract_splice_sites(gff_file, verbose=False):
    genes = dd(list)
    trans = {}

    for line in gff_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, phase, attributes = line.split('\t')
        except ValueError:
            continue

        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        attr_dict = parse_attributes_gff3(attributes)

        # Try to get transcript ID (Parent usually refers to mRNA/transcript)
        transcript_id = attr_dict.get('Parent')
        gene_id = attr_dict.get('gene_id') or attr_dict.get('gene') or transcript_id

        if transcript_id is None:
            continue  # skip if no parent/transcript ID

        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[gene_id].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Merge nearby exons (≤5bp apart)
    for tran, (chrom, strand, exons) in trans.items():
        exons.sort()
        merged = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - merged[-1][1] <= 5:
                merged[-1][1] = exons[i][1]
            else:
                merged.append(exons[i])
        trans[tran] = [chrom, strand, merged]

    # Print junctions
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, exons[i-1][1], exons[i][0], strand))

    for chrom, left, right, strand in sorted(junctions):
        print(f"{chrom}\t{left-1}\t{right-1}\t{strand}")  # zero-based

    # Verbose stats
    if verbose:
        exon_lengths, intron_lengths, trans_lengths = Counter(), Counter(), Counter()
        for chrom, strand, exons in trans.values():
            tran_len = 0
            for i, exon in enumerate(exons):
                exon_len = exon[1] - exon[0] + 1
                exon_lengths[exon_len] += 1
                tran_len += exon_len
                if i > 0:
                    intron_lengths[exon[0] - exons[i-1][1]] += 1
            trans_lengths[tran_len] += 1

        print(f"genes: {len(genes)}, genes with multiple isoforms: {sum(len(v)>1 for v in genes.values())}", file=stderr)
        print(f"transcripts: {len(trans)}, transcript avg. length: {sum(trans_lengths.elements())//len(trans)}", file=stderr)
        print(f"exons: {sum(exon_lengths.values())}, exon avg. length: {sum(exon_lengths.elements())//sum(exon_lengths.values())}", file=stderr)
        print(f"introns: {sum(intron_lengths.values())}, intron avg. length: {sum(intron_lengths.elements())//sum(intron_lengths.values())}", file=stderr)
        print(f"average number of exons per transcript: {sum(exon_lengths.values())//len(trans)}", file=stderr)

if __name__ == '__main__':
    parser = ArgumentParser(description="Extract splice sites from a GFF3 file (for HISAT2 --ss input)")
    parser.add_argument('gff_file', nargs='?', type=FileType('r'), help='Input GFF3 file (use "-" for stdin)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print statistics to stderr')

    args = parser.parse_args()
    if not args.gff_file:
        parser.print_help()
        exit(1)

    extract_splice_sites(args.gff_file, args.verbose)

