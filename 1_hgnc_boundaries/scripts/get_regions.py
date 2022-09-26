#!/usr/bin/env python3

RECORD_FIELDS = (
    'contig',
    'source',
    'type',
    'start',
    'end',
    'unknown_1',
    'strand',
    'unknown_2',
    'annotations',
)

def main():
    hgnc_ids = dict()
    with open('1_regions/hgnc_ids.txt', 'r') as fh:
        line_token_gen = (line.rstrip().split(' ') for line in fh)
        for symbol, hgnc_id in line_token_gen:
            hgnc_ids[hgnc_id] = symbol

    gff_records = list()
    with open('data/gencode.v41.chr_patch_hapl_scaff.annotation.gff3', 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line_tokens = line.split('\t')
            record = {k: v for k, v in zip(RECORD_FIELDS, line_tokens)}

            if record['type'] != 'exon':
                continue

            anns = dict()
            anns_tokens = record['annotations'].split(';')
            for ats in anns_tokens:
                k, v = ats.split('=')
                assert k not in anns_tokens
                anns[k] = v

            if 'hgnc_id' not in anns:
                continue


            if anns['hgnc_id'] in hgnc_ids:
                gene = hgnc_ids[anns['hgnc_id']].replace('_', '-').upper()
                anns_str = f'gene:{gene};type:exon;ensembl_id:{anns["ID"]}'
                data = (record['contig'], record['start'], record['end'], anns_str)
                print(*data, sep='\t')


if __name__ == '__main__':
    main()
