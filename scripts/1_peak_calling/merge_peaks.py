import os
import pandas as pd
import argparse
import numpy as np
from collections import defaultdict

CHROMS = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']
# Retrieved from SGD_R64-3-1_genome.fna
chr_lengths = {
                'chrI': 230218,
                'chrII': 813184,
                'chrIII': 316620,
                'chrIV': 1531933,
                'chrV': 576874,
                'chrVI': 270161,
                'chrVII': 1090940,
                'chrVIII': 562643,
                'chrIX': 439888,
                'chrX': 745751,
                'chrXI': 666816,
                'chrXII': 1078177,
                'chrXIII': 924431,
                'chrXIV': 784333,
                'chrXV': 1091291,
                'chrXVI': 948066,
                'chrMT': 85779
            }

def get_already_merged(output_dir: str) -> dict:
    '''Get groups which replicates are already merged.
    
    Returns:
        Dict{strain:[str(group1), str(group2), ..]}
    '''
    already_merged = defaultdict(list)

    for strain in os.listdir(output_dir):
        strain_path = os.path.join(output_dir, strain)

        if os.path.isdir(strain_path):
            already_merged[strain] = [sample.split('.')[0] for sample in os.listdir(strain_path)
                                        if sample.endswith('.bed')]
        else:
            continue
    return already_merged

def cut_telomers(peaks: pd.DataFrame, telomer_len: int) -> pd.DataFrame:
    '''Function to filter out peaks that are in the first and last `telomer_len` bp of each chromosome.

    Returns:
        Filtered pd.DataFrame
    '''
    filtered_peaks_list = [
        chr_peaks[(chr_peaks['start'] > telomer_len) & 
                  (chr_peaks['end'] < (chr_lengths[chr] - telomer_len))]
        for chr, chr_peaks in peaks.groupby('chr')
        if chr in chr_lengths
    ]
    
    # Filter out empty DataFrames and concatenate
    filtered_peaks_list = [df for df in filtered_peaks_list if not df.empty]
    
    if filtered_peaks_list:
        return pd.concat(filtered_peaks_list, ignore_index=True)
    return pd.DataFrame(columns=peaks.columns)

def load_replicates(paths_list: list[str], filter_telomers: bool, telomer_len: int) -> dict:
    '''Load .bed to pd.DataFrames for corresponding replicate from "paths_list"
    If specified, filter out peaks inside telomer region.

    Returns:
        Dict{'A':pd.DataFrame,'B':pd.DataFrame,'C':pd.DataFrame}
    '''
    replicates = defaultdict(pd.DataFrame)

    for rep_path in paths_list:
        if rep_path.endswith('.bed'):
            rep_filename = os.path.basename(rep_path)
            rep_name = rep_filename.split('.')[0]
            rep = rep_name.split('_')[-1]

            if rep in ['A', 'B', 'C']:
                # Read the BED file and ensure numeric types for positions
                df = pd.read_csv(rep_path, usecols=[0,1,2,3], sep='\t', header=None)
                df.columns = ['chr', 'start', 'end', 'center']
                
                # Convert numeric columns to appropriate types
                df['start'] = df['start'].astype(int)
                df['end'] = df['end'].astype(int)
                df['center'] = df['center'].astype(float)
                
                # Cut telomers if specified
                if filter_telomers:
                    df = cut_telomers(df, telomer_len)
                    
                replicates[rep] = df
            else:
                print("Warning: replicate ID is not 'A', 'B', 'C'")
                return False
    return replicates

def two_out_of_three(replicates: dict) -> pd.DataFrame:
    '''Algorithm to find overlapping peaks in at least 2 out of 3 replicates.
    Works for 3 samples only.

    Returns:
        pd.DataFrame
    '''
    merged = {'chr': [], 'start': [], 'end': [], 'center': []}

    for chrom in CHROMS:
        # Collect all intervals with their replicate ID
        intervals = []
        for rep, rep_df in replicates.items():
            rep_chrom = rep_df[rep_df['chr'] == chrom]

            for _, row in rep_chrom.iterrows():
                intervals.append((row['start'], row['end'], row['center'], rep))
        
        if not intervals:
            print(f'Warning: chromosome {chrom} has no peaks for all replicates.')
            continue

        # Merging algorithm    
        # Sort by start position
        intervals.sort(key=lambda x: x[0])

        current_start = None
        current_end = None
        active_reps = set()

        for start, end, center, rep in intervals:
            # Init start and end
            if current_start is None and current_end is None:
                current_start = start
                current_end = end
                active_reps.add((center, rep))

            # Overlap, new start inside, new end outside
            elif start < current_end and end >= current_end: 
                # Change only start
                current_start = start
                active_reps.add((center, rep))

            # Overlap, new interval is fully inside 'current' interval
            elif start >= current_end and end <= current_end: 
                current_start = start
                current_end = end
                active_reps.add((center,rep))
            
            # Overlap, start outside, end inside
            elif start <= current_start and end < current_end: 
                # Change only end
                current_end = end
                active_reps.add((center,rep))

            elif start > current_end: #leaving overlapping
                # We want record only overlapped regions
                if len(active_reps) >= 2 and current_start is not None:
                    center_mean = int(np.mean([rep[0] for rep in active_reps]))
                    # Found an overlap region
                    # Record new overlap region using mean of centers of overlapped intervals
                    # Start is -75 bp from new center and end is +75 bp.
                    merged['chr'].append(chrom)
                    merged['start'].append(int(center_mean - 75))
                    merged['end'].append(int(center_mean + 75))
                    merged['center'].append(center_mean)

                # Reset for the new interval
                current_start = start
                current_end = end
                active_reps = {(center, rep)}
        
        # Check for the last interval
        if len(active_reps) >= 2 and current_start is not None:
            center_mean = np.mean([rep[0] for rep in active_reps])
            merged['chr'].append(chrom)
            merged['start'].append(int(center_mean - 75))
            merged['end'].append(int(center_mean + 75))
            merged['center'].append(center_mean)
    
    return pd.DataFrame(merged)

def save_merged(merged: pd.DataFrame, output_dir: str, group: str):
    '''Save merged df in .bed file to provided output directory
    '''
    output_bed = os.path.join(output_dir, group + '.bed')
    merged.to_csv(output_bed, sep='\t', index=False, header=False)

    # Create summary .csv of peaks count for different experiments for given strain
    peaks_count = merged.shape[0]
    peaks_count_row = pd.DataFrame(
        {'sample': [group],
        'peaks': [peaks_count]})
    peaks_count_row.to_csv(os.path.join(output_dir, 'peaks_number.csv'), mode='a', index=False, header=not os.path.exists(os.path.join(output_dir, 'peaks_number.csv')))


########################### PARSE ARGUMENTS AND RUN MERGING ########################### 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str,
                    help="Path to directory of directories corresponding to each strain with .bed files from different experiments and all three replicates.")
    parser.add_argument("-o", "--output_dir", type=str,
                    help="Path to directory of directories corresponding to each strain with .bed files merged across replicates from one experiment.")
    parser.add_argument("-t", "--telomer_length", type=int, default=1000,
                    help="Provide telomer length to filter first and last bp for each chromosome. Default 1000 bp.")
    parser.add_argument("-f", "--filter_telomers", type=bool, default=True,
                    help="If applied, filter out telomers with specified length. Default True.")
    args = parser.parse_args()

    # Check arguments for telomers
    if not args.filter_telomers and args.telomer_length:
        print("Warning: Cut telomers set as False, telomer length is ignored.")
    else:
        print(f"Cut telomers set as True (default), peaks in first and last {args.telomer_length} bp of each chromosome are ignored.\n")
        print("If needed, change this by setting -f {True/False} -t {number of bp}.")

    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # If some samples already merged in output directory, skip them
    to_exclude = get_already_merged(args.output_dir)
    to_merge = defaultdict(lambda: defaultdict(list))

    for strain in os.listdir(args.input_dir):
        strain_input_dir = os.path.join(args.input_dir,strain)
        for rep in os.listdir(strain_input_dir):
            if rep.endswith('.bed'):
                group = '_'.join(rep.split('_')[:-1]) # get filename without replicate number - this is grouping name
                if group in to_exclude[strain]:
                    print(f'Strain {strain} with group name {group} is already merged')
                else:
                    to_merge[strain][group].append(os.path.join(strain_input_dir, rep))
    
    # Merge replicates within defined group
    for strain in to_merge:
        if not os.path.exists(os.path.join(args.output_dir, strain)):
            os.makedirs(os.path.join(args.output_dir, strain))

        for group in to_merge[strain]:
            # Load replicates .bed file and merge them
            if len(to_merge[strain][group]) == 3:
                replicates = load_replicates(to_merge[strain][group],
                                             filter_telomers=args.filter_telomers,
                                             telomer_len=args.telomer_length)
                if replicates:
                    merged = two_out_of_three(replicates)
                    save_merged(merged, os.path.join(args.output_dir, strain), group)
                    print(f'Merged peaks for {strain} {group}')
                else:
                    print(f"Warning: Check replicates filenames for {group}")
            else:
                print(f"Warning: Incomplete set of replicates for {group}")
                continue

if __name__ == "__main__":
    main()