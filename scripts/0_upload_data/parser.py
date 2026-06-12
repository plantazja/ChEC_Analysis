import re
import pandas as pd
import argparse
import os
from os import listdir
from os.path import isfile, join, isdir

def parse_per_row(row_raw, df_new):
    pool, sample, sample_name = row_raw.split(',')

    # Exctract merge group
    merge_match = re.search(r"^(.*?)_\w{1}$", sample_name)
    if not merge_match:
        raise ValueError(f"Warning: Could not extract merge group from {sample_name}")
    else:
        merge_group = merge_match.group(1)
    
    archive_path = "/archive/donczew/HighThroughputSequencing/"
    path_to_pool = os.path.join(archive_path, pool)
    
    # Check if directory exist
    if not os.path.isdir(path_to_pool):
        raise ValueError(f"No directory {path_to_pool}")
    
    dir_after_pool = [d for d in listdir(path_to_pool) if isdir(join(path_to_pool, d))]
    if not dir_after_pool:
        raise ValueError(f"No subdirectories found in {path_to_pool}")
    
    path_to_fastqs = os.path.join(path_to_pool, dir_after_pool[0])

    # Check if file exist
    onlyfiles = [f for f in listdir(path_to_fastqs) if isfile(join(path_to_fastqs, f))]
    onlyfastqs = [f for f in onlyfiles if f.endswith('.fastq.gz')]

    fastq1_fnam = None
    fastq2_fnam = None

    for fnam in onlyfastqs:
        if sample == fnam.split('_')[0]:
            if "R1" == fnam.split('_')[2]:
                fastq1_fnam = os.path.join(path_to_fastqs, fnam)

            elif "R2" == fnam.split('_')[2]:
                fastq2_fnam = os.path.join(path_to_fastqs, fnam)
        
    if not fastq1_fnam or not fastq2_fnam:
        raise ValueError(f"Missing FASTQ files for sample {sample} in {path_to_fastqs}")
    
    new_row = {
        'sample': sample_name,
        'fastq1': fastq1_fnam,
        'fastq2': fastq2_fnam,
        'merge_group': merge_group
    }
    df_new = pd.concat([df_new, pd.DataFrame([new_row])], ignore_index=True)
    return df_new

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="takes as input path to csv file with columns: 'pool', 'sample', 'sample_name")
    parser.add_argument('-o', '--output_dir', help='takes as output path to outpu directory')
    args = parser.parse_args()
    
    df = pd.read_csv(args.input)
    try:
        list_df = (df['pool'] + ',' + df['sample'] + ',' + df['sample_name']).to_list()
    except Exception as e:
        print(f"{e}")

    df_new = pd.DataFrame(columns = ['sample', 'fastq1', 'fastq2', 'merge_group'])

    for row in list_df:
        try:
            df_new = parse_per_row(row, df_new)
        except Exception as e:
            print(f"{e}")

    pools = set(df['pool'])
    pools_num = [re.findall(r'\d+', p)[0] for p in pools]
    pools_num_str = '_'.join(pools_num)
    output_fnam = 'RDHTSP' + pools_num_str + '.csv'
    output_path = os.path.join(args.output_dir, output_fnam) 

    df_new.to_csv(output_path, index=False)

if __name__ == '__main__':
    main()