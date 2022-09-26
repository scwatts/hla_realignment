#!/usr/bin/env python3
import csv


def main():
    with open('data/hg38.chromAlias.txt', 'r') as fh:
        line_token_gen = csv.reader(fh, delimiter='\t')

        header_tokens = next(line_token_gen)
        header_tokens[0] = header_tokens[0].replace('# ', '')

        records = dict()
        for line_tokens in line_token_gen:
            record = {k: v for k, v in zip(header_tokens, line_tokens)}
            record['genbank'] not in records

            assert record['genbank'] not in records
            records[record['genbank']] = record

    with open('1_regions/regions_gencode_unmerged.bed', 'r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for line_tokens in line_token_gen:

            if line_tokens[0] == 'chr6':
                print(*line_tokens, sep='\t')
                continue

            assert line_tokens[0] in records

            record = records[line_tokens[0]]
            line_tokens[0] = record['ucsc']
            print(*line_tokens, sep='\t')


if __name__ == '__main__':
    main()
