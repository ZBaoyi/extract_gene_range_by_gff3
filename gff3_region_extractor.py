#!/usr/bin/env python3
import argparse

def parse_attributes(attr_str):
    """
    Parse the GFF3 attribute string and return a dictionary.
    """
    attributes = {}
    for part in attr_str.strip().split(';'):
        if part:
            kv = part.split('=', 1)
            if len(kv) == 2:
                key, value = kv
                attributes[key.strip()] = value.strip()
    return attributes

def merge_intervals(intervals):
    """
    Merge overlapping or contiguous intervals.
    Returns a list of merged intervals.
    """
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1] + 1:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def calculate_relative_position(start, end, gene_record):
    """
    Calculate the relative position based on the strand.
    For positive strand, relative position is just the position.
    For negative strand, relative position is the reverse.
    """
    if gene_record['strand'] == '+':
        return start - gene_record['start'] + 1, end - gene_record['start'] + 1
    else:
        return gene_record['end'] - start + 1, gene_record['end'] - end + 1


def main():
    parser = argparse.ArgumentParser(
        description="Extract specified region genomic coordinates for a gene from a GFF3 file.\npython gff3_region_extractor.py -g LOC_Os01g01010 -f path/to/your.gff3 -r promoter -p 2000 -o output.txt"
    )
    parser.add_argument('-g', '--gene', required=True,
                        help='Gene ID or name, e.g., LOC_Os01g01010')
    parser.add_argument('-f', '--file', required=True,
                        help='Path to the GFF3 file')
    parser.add_argument('-r', '--region', required=True,
                        choices=['promoter', 'exon', 'intron', 'five_prime_UTR', 'three_prime_UTR', 'CDS', 'gene_body'],
                        help='Region type to extract: promoter, exon, intron, five_prime_UTR, three_prime_UTR, CDS, or gene_body')
    parser.add_argument('-p', '--promoter', type=int, default=2000,
                        help='Promoter length when region is promoter (default: 2000)')
    parser.add_argument('-o', '--output',
                        help='Output file to save results. If not provided, results will be printed to stdout.')
    args = parser.parse_args()

    gene_query = args.gene
    gff3_file = args.file
    region_type = args.region
    promoter_length = args.promoter
    output_file = args.output

    gene_record = None
    mrna_ids = []  # List to store mRNA IDs associated with the gene
    # Dictionary to store features (lists of intervals)
    features = {
        'exon': [],
        'CDS': [],
        'five_prime_UTR': [],
        'three_prime_UTR': []
    }

    output_lines = []

    # Parse the GFF3 file
    with open(gff3_file) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_str = parts
            start = int(start)
            end = int(end)
            attrs = parse_attributes(attr_str)

            if ftype == "gene":
                gene_id = attrs.get("ID", "")
                gene_name_attr = attrs.get("Name", "")
                # Match gene by ID or Name
                if gene_id == gene_query or gene_name_attr == gene_query:
                    gene_record = {
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'attributes': attrs
                    }
            elif ftype == "mRNA":
                parent = attrs.get("Parent", "")
                if gene_record is not None:
                    if parent == gene_record['attributes'].get("ID", ""):
                        mrna_ids.append(attrs.get("ID"))
                else:
                    if gene_query in parent:
                        mrna_ids.append(attrs.get("ID"))
            else:
                # Collect features if they belong to the relevant mRNA(s)
                parent = attrs.get("Parent", "")
                parent_list = parent.split(',') if parent else []
                if any(mid in parent_list for mid in mrna_ids):
                    if ftype in features:
                        features[ftype].append((start, end))

    if gene_record is None:
        msg = f"Gene {gene_query} not found in {gff3_file}."
        if output_file:
            with open(output_file, 'w') as fout:
                fout.write(msg + "\n")
        else:
            print(msg)
        return

    if region_type == "gene_body":
        rel_start, rel_end = calculate_relative_position(gene_record['start'], gene_record['end'], gene_record)
        output_lines.append(f"{gene_record['seqid']}\t{gene_record['start']}\t{gene_record['end']}\t{gene_query}\t{rel_start}-{rel_end}")
    elif region_type == "promoter":
        if gene_record['strand'] == '+':
            promoter_start = max(1, gene_record['start'] - promoter_length)
            promoter_end = gene_record['start'] - 1
        else:
            promoter_start = gene_record['end'] + 1
            promoter_end = gene_record['end'] + promoter_length
        rel_start, rel_end = calculate_relative_position(promoter_start, promoter_end, gene_record)
        output_lines.append(f"{gene_record['seqid']}\t{promoter_start}\t{promoter_end}\t{gene_query}-Promoter({promoter_length})\t{gene_record['strand']}\t{rel_start}-{rel_end}")
    elif region_type == "exon":
        exons = merge_intervals(features['exon'])
        if not exons:
            output_lines.append("# No exon records found.")
        else:
            if gene_record['strand'] == '+':
                for i, exon in enumerate(exons):
                    rel_start, rel_end = calculate_relative_position(exon[0], exon[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{exon[0]}\t{exon[1]}\t{gene_query}-Exon-{i+1}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
            else:
                num_exons = len(exons)
                for i, exon in enumerate(exons):
                    rel_start, rel_end = calculate_relative_position(exon[0], exon[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{exon[0]}\t{exon[1]}\t{gene_query}-Exon-{num_exons-i}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
    elif region_type == "intron":
        exons = merge_intervals(features['exon'])
        if len(exons) < 2:
            output_lines.append("# Not enough exon information to calculate introns.")
        else:
            introns = []
            for i in range(len(exons) - 1):
                intron_start = exons[i][1] + 1
                intron_end = exons[i+1][0] - 1
                if intron_start <= intron_end:
                    introns.append((intron_start, intron_end))
            if not introns:
                output_lines.append("# No intron found.")
            else:
                if gene_record['strand'] == '+':
                    for i, intron in enumerate(introns):
                        rel_start, rel_end = calculate_relative_position(intron[0], intron[1], gene_record)
                        output_lines.append(f"{gene_record['seqid']}\t{intron[0]}\t{intron[1]}\t{gene_query}-Intron-{i+1}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
                else:
                    num_introns = len(introns)
                    for i, intron in enumerate(introns):
                        rel_start, rel_end = calculate_relative_position(intron[0], intron[1], gene_record)
                        output_lines.append(f"{gene_record['seqid']}\t{intron[0]}\t{intron[1]}\t{gene_query}-Intron-{num_introns-i}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
    elif region_type == "five_prime_UTR":
        five_utrs = sorted(features['five_prime_UTR'], key=lambda x: x[0])
        if not five_utrs:
            output_lines.append("# No 5' UTR records found.")
        else:
            if gene_record['strand'] == '+':
                for i, utr in enumerate(five_utrs):
                    rel_start, rel_end = calculate_relative_position(utr[0], utr[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{utr[0]}\t{utr[1]}\t{gene_query}-5UTR-{i+1}\t{gene_record['strand']}{rel_start}-{rel_end}")
            else:
                num_utrs = len(five_utrs)
                for i, utr in enumerate(five_utrs):
                    rel_start, rel_end = calculate_relative_position(utr[0], utr[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{utr[0]}\t{utr[1]}\t{gene_query}-5UTR-{num_utrs-i}\t{gene_record['strand']}{rel_start}-{rel_end}")
    elif region_type == "three_prime_UTR":
        three_utrs = sorted(features['three_prime_UTR'], key=lambda x: x[0])
        if not three_utrs:
            output_lines.append("# No 3' UTR records found.")
        else:
            if gene_record['strand'] == '+':
                for i, utr in enumerate(three_utrs):
                    rel_start, rel_end = calculate_relative_position(utr[0], utr[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{utr[0]}\t{utr[1]}\t{gene_query}-3UTR-{i+1}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
            else:
                num_utrs = len(three_utrs)
                for i, utr in enumerate(three_utrs):
                    rel_start, rel_end = calculate_relative_position(utr[0], utr[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{utr[0]}\t{utr[1]}\t{gene_query}-3UTR-{num_utrs-i}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
    elif region_type == "CDS":
        cds_list = sorted(features['CDS'], key=lambda x: x[0])
        if not cds_list:
            output_lines.append("# No CDS records found.")
        else:
            if gene_record['strand'] == '+':
                for i, cds in enumerate(cds_list):
                    rel_start, rel_end = calculate_relative_position(cds[0], cds[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{cds[0]}\t{cds[1]}\t{gene_query}-CDS-{i+1}\t{gene_record['strand']}\t{rel_start}-{rel_end}")
            else:
                num_cdss = len(cds_list)
                for i, cds in enumerate(cds_list):
                    rel_start, rel_end = calculate_relative_position(cds[0], cds[1], gene_record)
                    output_lines.append(f"{gene_record['seqid']}\t{cds[0]}\t{cds[1]}\t{gene_query}-CDS-{num_cdss-i}\t{gene_record['strand']}\t{rel_start}-{rel_end}")

    result = "\n".join(output_lines)

    if output_file:
        with open(output_file, 'w') as fout:
            fout.write(result + "\n")
        print(f"Results saved to {output_file}")
        print("-"*30)
        print(result)
        print("-"*30)
    else:
        print(result)

if __name__ == "__main__":
    main()
