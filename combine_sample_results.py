import sys
import argparse
import pandas as pd
import os

def combine(samples):
    df = {}
    for sample in samples:
        s = pd.Series.from_csv(sample, sep='\t', header=None).to_dict()
        df[os.path.basename(sample)] = s
    df = pd.DataFrame(df)
    df = df.fillna(value=0)
    for c in df.columns:
        df[c] = df[c].astype(int)
    return df

if __name__ == '__main__':
    # parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_1", help="sample 1 file",type=str)
    parser.add_argument("other_samples", help="list of the other samples files to combine", nargs='+', type=str)
    parser.add_argument("-o", "--output_file", help="the combined output file", type=str, default='combined_samples.txt')
    args = parser.parse_args()
    samples = [args.sample_1] + [x for x in args.other_samples]
    df = combine(samples)
    df.to_csv(args.output_file, sep='\t')
