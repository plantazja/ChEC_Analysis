import sys, os
import os.path
import pandas as pd
import numpy as np

def process_bed(in_path, out_path):
    '''
    Parsing raw .bed output from HOMER to new .bed file with 'center' column.
    '''
    chroms = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']

    bed_df = pd.read_csv(in_path, sep='\t', header=None, usecols=[1,2,3,11], names=['chr','start','end','p-value vs Control'], comment='#')
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
    peak_size = bed_df['end'] - bed_df['start']
    bed_df['center'] = bed_df['start'] + peak_size / 2
    bed_df['center'] = bed_df['center'].astype(int)
    bed_df.chr = pd.Categorical(values=bed_df.chr, categories=chroms, ordered=True)
    bed_df.sort_values(['chr', 'center'], inplace=True)

    # Remove peaks that mapped on rDNA locus on chrXII
    bed_df = bed_df.drop(bed_df[(bed_df['chr'] == 'chrXII') &
                                (bed_df['end'] >= 451000) &
                                (bed_df['start'] <= 469000)].index)
    
    bed_df[['chr', 'start', 'end', 'center', 'p-value vs Control']].to_csv(out_path, header=False, index=False, sep='\t')

    return bed_df.shape[0]  # Return number of peaks for this sample

def main():
    beds_dir = sys.argv[1]
    out_dir = sys.argv[2] 

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Init dictionary of pd.DataFrames per strain.
    # each row in df corresponds to separate sample.
    dict_peaks = {}
    
    for f in os.listdir(beds_dir):
        if f.endswith('.bed'):
            in_path = os.path.join(beds_dir, f)

            # Exctract all covariants
            sample_name = f.replace('.bed', '')
            strain = sample_name.split('_')[0]

            # Get name for new .bed            
            processed_bed = f
            
            # Create strain directory in output directory
            out_path = os.path.join(out_dir, strain)
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            out_path = os.path.join(out_path, processed_bed)

            # Skip file if already processed
            if os.path.isfile(out_path):
                print(f'File for {sample_name} is already exists in {out_path}, skipping..')
                continue
            else:            
                # Get peak count per sample and save info in row
                peaks_count = process_bed(in_path, out_path)
                new_row = pd.DataFrame({'sample': [sample_name], 'peaks': [peaks_count]})
                
                # Append row to corresponding strain pd.DataFrame
                if strain not in dict_peaks:
                    dict_peaks[strain] = pd.DataFrame(columns=['sample', 'peaks'])
                    dict_peaks[strain] = pd.concat([dict_peaks[strain], new_row], ignore_index=True)
                else:
                    dict_peaks[strain] = pd.concat([dict_peaks[strain], new_row], ignore_index=True)

    for strain, df_peaks in dict_peaks.items():
        #Sort the DataFrame by sample name, just to have some order in .csv file        
        df_peaks.sort_values(by='sample', inplace=True)

        # Save peaks_number.csv for all samples for given strain in strain directory 
        out_path = os.path.join(out_dir, strain)
        out_path = os.path.join(out_path, 'peaks_number.csv')
        df_peaks.to_csv(out_path, mode='a', index=False, header=not os.path.exists(out_path))

if __name__ == "__main__":
    main()