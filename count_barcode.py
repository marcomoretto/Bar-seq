import sys
import itertools
import Levenshtein
import argparse

def DEPRECATED__get_barcode_to_gene_mapping(file, gene_column, barcode_column, unmatched_file):
    barcode_to_gene_mapping = {}
    header = []
    unmatched = []
    with open(file) as f:
        header = [x.strip() for x in f.next().strip().split('\t')]
        orf_idx = header.index(gene_column)
        barcode_idx = header.index(barcode_column)
        f.next()
        for l in f:
            try:
                s = [x.strip() for x in l.strip().split('\t')]
                orf = s[orf_idx]
                barcode = s[barcode_idx]
                if barcode in barcode_to_gene_mapping and orf != barcode_to_gene_mapping[barcode]:
                    sys.exit('\n\n*********\n*ERROR!\n*********\n\n' + barcode + ' is mapping to more than one gene: ' + orf + ' and ' + barcode_to_gene_mapping[barcode] + '\n\nABORTED!\n')
                barcode_to_gene_mapping[barcode] = orf
            except Exception as e:
                unmatched.append(l.strip())
                pass
    with open(unmatched_file, 'w') as f:
        f.write('\t'.join(header))
        f.write('\n')
        f.write('\n'.join(unmatched))
    return barcode_to_gene_mapping

def get_barcode_to_gene_mapping(file, gene_column, barcode_column, barcode_up_down, unmatched_file):
    barcode_to_gene_mapping = {}
    header = []
    unmatched = []
    with open(file) as f:
        header = [x.strip() for x in f.next().strip().split('\t')]
        orf_idx = header.index(gene_column)
        barcode_idx = header.index(barcode_column)
        for l in f:
            try:
                s = [x.strip() for x in l.strip().split('\t')]
                if s[0][-1:] == barcode_up_down:
                    orf = s[orf_idx]
                    barcode = s[barcode_idx]
                    if barcode in barcode_to_gene_mapping and orf != barcode_to_gene_mapping[barcode]:
                        sys.exit('\n\n*********\n*ERROR!\n*********\n\n' + barcode + ' is mapping to more than one gene: ' + orf + ' and ' + barcode_to_gene_mapping[barcode] + '\n\nABORTED!\n')
                    barcode_to_gene_mapping[barcode] = orf
            except Exception as e:
                unmatched.append(l.strip())
                pass
    with open(unmatched_file, 'w') as f:
        f.write('\t'.join(header))
        f.write('\n')
        f.write('\n'.join(unmatched))
    return barcode_to_gene_mapping

def get_total_lines(file):
    with open(file) as f:
        num_lines = sum(1 for line in f)
    return num_lines

def print_results(output_file, barcode_to_gene_mapping, barcode_counts, output_unmatched):
    known_barcode = set.intersection(set(barcode_to_gene_mapping.keys()), set(barcode_counts.keys()))
    unknown_barcode = set(barcode_counts.keys()) - known_barcode
    tot_counts = {k:0 for k in barcode_to_gene_mapping.values()}
    for barcode,gene in barcode_to_gene_mapping.iteritems():
        tot_counts[gene] += barcode_counts[barcode]
    with open(output_file, 'w') as f:
        for gene, counts in tot_counts.iteritems():
            f.write(gene + '\t' + str(counts) + '\n')
        if output_unmatched:
            for barcode in unknown_barcode:
                f.write(barcode + '\t' + str(barcode_counts[barcode]) + '\n')

def count_barcode(input_file, index_len, primer_len, barcode_len, total_seq, barcode_to_gene_mapping):
    barcode_counts = {k:0 for k in barcode_to_gene_mapping.keys()}
    done = 0.0
    with open(input_file) as f:
        for line in itertools.islice(f, 1, None, 4):
            done +=1
            barcode = line.strip()[index_len + primer_len:index_len + primer_len + barcode_len]
            if barcode not in barcode_counts:
                barcode_counts[barcode] = 1
            else:
                barcode_counts[barcode] += 1
            if done % 1000 == 0:
                msg = 'Counting ... ' + str(round(done / total_seq * 100, 2)) + '%'
                sys.stderr.write(msg)
                sys.stderr.flush()
                sys.stderr.write("\b" * (len(msg) + 1))
        sys.stderr.write("Counting ... done!      \n")
    return barcode_counts

def aggregate_counts(barcode_to_gene_mapping, barcode_counts, levenshtein_distance, multi_match_file):
    to_sum = {} # unknown_barcode -> gene_barcode
    multi_match = set()
    done = 0.0
    unknown_barcodes = set(barcode_counts.keys()) - set(barcode_to_gene_mapping.keys())
    for unknown_barcode in unknown_barcodes:
        done += 1
        for gene_barcode in barcode_to_gene_mapping.keys():
            if Levenshtein.distance(unknown_barcode, gene_barcode) <= levenshtein_distance:
                # check if the unknown_barcode is already associated with a gene
                # if that is the case we have to discard the unknown_barcode (it cannot map to more than one gene)
                if unknown_barcode in to_sum:
                    multi_match.add(unknown_barcode)
                    break
                else:
                    to_sum[unknown_barcode] = gene_barcode
        if done % 100 == 0:
            msg = 'Aggregating ... ' + str(round(done / float(len(unknown_barcodes)) * 100, 2)) + '%'
            sys.stderr.write(msg)
            sys.stderr.flush()
            sys.stderr.write("\b" * (len(msg) + 1))
    sys.stderr.write("Aggregating ... done!       \n")
    for barcode in multi_match:
        del barcode_counts[barcode]
    for unknown_barcode, gene_barcode in to_sum.iteritems():
        if unknown_barcode not in multi_match:
            barcode_counts[gene_barcode] += barcode_counts[unknown_barcode]
    for unknown_barcode in to_sum.keys():
        if unknown_barcode in barcode_counts:
            del barcode_counts[unknown_barcode]
    with open(multi_match_file, 'w') as f:
        f.write('\n'.join(list(multi_match)))

if __name__ == '__main__':
    # parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_file", help="fastq file",type=str)
    parser.add_argument("gene_barcode_mapping_file", help="file containing gene to barcode mapping",type=str)
    parser.add_argument("output_file", help="the count output file",type=str)
    parser.add_argument("-gc", "--gene_column_name", help="the column name for gene in gene_barcode_mapping_file (default Name),", type=str, default="Name")
    parser.add_argument("-bs", "--barcode_column_sequence", help="the column name for the barcode sequence in gene_barcode_mapping_file (default Sequence),", type=str, default="Sequence")
    parser.add_argument("-ud", "--barcode_up_down", help="use U or D (up or down) barcode in gene_barcode_mapping_file (default D),", type=str, default="D")
    #parser.add_argument("-bc", "--barcode_column_name", help="the column name for barcode in gene_barcode_mapping_file (default DNTAG_sequence_20mer),", type=str, default="DNTAG_sequence_20mer")
    parser.add_argument("-i", "--index_len", help="the index sequence len (default 5),", type=int, default=5)
    parser.add_argument("-p", "--primer_len", help="the primer sequence len (default 15),", type=int, default=15)
    parser.add_argument("-b", "--barcode_len", help="the barcode sequence len (default 20)", type=int, default=20)
    parser.add_argument("-l", "--levenshtein_distance", help="levenshtein distance (default 0)", type=int, default=0)
    parser.add_argument("-u", "--unmatched", help="file name with unmatched barcode to gene mapping (default unmatched_barcode.txt)", type=str, default='unmatched_barcode.txt')
    parser.add_argument("-ou", "--output_unmatched", help="print unmatched barcode in output_file (default FALSE)", type=bool, default=False)
    parser.add_argument("-m", "--multi_match", help="file name with multi-match barcodes to gene mapping (default multi_matched_barcode.txt)", type=str, default='multi_matched_barcode.txt')
    args = parser.parse_args()
    # get total number of sequences
    sys.stderr.write('Pre-processing ... ')
    tot_seq = float(get_total_lines(args.fastq_file)) / 4.0
    # read the barcode to gene file
    #barcode_to_gene_mapping = DEPRECATED__get_barcode_to_gene_mapping(args.gene_barcode_mapping_file, \
    #    args.gene_column_name, args.barcode_column_name, args.unmatched)
    barcode_to_gene_mapping = get_barcode_to_gene_mapping(args.gene_barcode_mapping_file, \
        args.gene_column_name, args.barcode_column_sequence, args.barcode_up_down, args.unmatched)
    sys.stderr.write('done!\n')
    # do the actual barcode count
    barcode_counts = count_barcode(args.fastq_file, args.index_len, args.primer_len, args.barcode_len, \
        tot_seq, barcode_to_gene_mapping)
    # check for multi-match and edit distance
    if args.levenshtein_distance > 0:
        aggregate_counts(barcode_to_gene_mapping, barcode_counts, args.levenshtein_distance, args.multi_match)
    # write results to file
    print_results(args.output_file, barcode_to_gene_mapping, barcode_counts, args.output_unmatched)
    sys.stderr.write('All done! \n')
