usage: gff3_region_extractor.py [-h] -g GENE -f FILE -r
                                {promoter,exon,intron,five_prime_UTR,three_prime_UTR,CDS,gene_body}
                                [-p PROMOTER] [-o OUTPUT]

Extract specified region genomic coordinates for a gene from a GFF3 file.
python gff3_region_extractor.py -g LOC_Os01g01010 -f path/to/your.gff3 -r
promoter -p 2000 -o output.txt

options:
  -h, --help            show this help message and exit
  -g GENE, --gene GENE  Gene ID or name, e.g., LOC_Os01g01010
  -f FILE, --file FILE  Path to the GFF3 file
  -r {promoter,exon,intron,five_prime_UTR,three_prime_UTR,CDS,gene_body}, --region {promoter,exon,intron,five_prime_UTR,three_prime_UTR,CDS,gene_body}
                        Region type to extract: promoter, exon, intron,
                        five_prime_UTR, three_prime_UTR, CDS, or gene_body
  -p PROMOTER, --promoter PROMOTER
                        Promoter length when region is promoter (default:
                        2000)
  -o OUTPUT, --output OUTPUT
                        Output file to save results. If not provided, results
                        will be printed to stdout.
