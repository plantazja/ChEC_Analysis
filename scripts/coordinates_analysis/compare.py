import pandas as pd
import numpy as np

CHROMS = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']

def compare_samples(samples_dict: dict, include_neighbours=0) -> pd.DataFrame:
    samples_name = [str(s) for s in samples_dict.keys()]

    # Init columns names including columns with samples names
    cols = ['chr', 'start', 'end', 'center'] + samples_name
    
    # Init dict for tracking overlaps with other samples
    overlap_tracker = {col: [] for col in cols}

    # Collect all intervals and then sort them in ascending order
    for chr in CHROMS:
        intervals = []
        for sample, s_df in samples_dict.items():
            s_chr = s_df[s_df['chr'] == chr]
            for _, row in s_chr.iterrows():
                intervals.append((row['start'], row['end'], row['center'], sample))

        if not intervals:  # Skip empty chromosomes
            continue

        # Sort by start position
        intervals.sort(key=lambda x: x[0])

        # Init current start and end for overlap region
        current_start = intervals[0][0]
        current_end = intervals[0][1]
        current_centers = [intervals[0][2]]  # Store centers as list
        active_samples = {intervals[0][3]}   # Store just sample names as set

        for start, end, center, sample in intervals[1:]:
            if current_start <= (end + include_neighbours) and current_end >= (start - include_neighbours):  # Overlap!
                current_start = max(current_start, start)
                current_end = min(current_end, end)
                current_centers.append(center)
                active_samples.add(sample)
            else:  # No overlap: record previous region
                # Calculate mean center (rounded)
                center_mean = int(round(np.mean(current_centers)))

                overlap_tracker['chr'].append(chr)
                overlap_tracker['start'].append(int(center_mean - 75))
                overlap_tracker['end'].append(int(center_mean + 75))
                overlap_tracker['center'].append(center_mean)
                
                # Mark which samples are in this region
                for s in samples_name:
                    if s in active_samples:
                        overlap_tracker[s].append(1)
                    else:
                        overlap_tracker[s].append(0)
                
                # Reset for the new interval
                current_start = start
                current_end = end
                current_centers = [center]
                active_samples = {sample}
        
        # Add last region
        center_mean = int(round(np.mean(current_centers)))
        overlap_tracker['chr'].append(chr)
        overlap_tracker['start'].append(int(center_mean - 75))
        overlap_tracker['end'].append(int(center_mean + 75))
        overlap_tracker['center'].append(center_mean)
        
        for s in samples_name:
            if s in active_samples:
                overlap_tracker[s].append(1)
            else:
                overlap_tracker[s].append(0)

    return pd.DataFrame(overlap_tracker)

def compare_pairwise(origins_dict, factor_dict, include_neighbours=0):
    ''' Input as dictionary {'ARS': ars_df}, {'factor_name': factor_df} '''
    
    CHROMS = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']
    
    origins_name = list(origins_dict.keys())[0]
    origin_df = origins_dict[origins_name]

    factor_name = list(factor_dict.keys())[0]
    factor_df = factor_dict[factor_name]

    # create tracker with same index row as origin_df
    overlap_tracker = pd.DataFrame(index=origin_df.index, columns=[f'{factor_name}_start',
                                                                    f'{factor_name}_end',
                                                                    f'{factor_name}_center'])

    for i, origin_row in origin_df.iterrows():
        curr_overlap = factor_df[((factor_df['end'] + include_neighbours) >= origin_row['start']) & # extend length of factor to find neighboring origins 
                                ((factor_df['start'] - include_neighbours) <= origin_row['end']) & # extend length of factor to find neighboring origins 
                                (factor_df['chr'] == origin_row['chr'])]
        # If no overlap for current row
        if curr_overlap.empty:
            overlap_tracker.loc[i, f'{factor_name}_start'] = np.nan
            overlap_tracker.loc[i, f'{factor_name}_end'] = np.nan
            overlap_tracker.loc[i, f'{factor_name}_center'] = np.nan
        elif len(curr_overlap) == 1: # if only one row, meaning only one interval is overlapping with ARS
            overlap_row = curr_overlap.iloc[0]

            overlap_tracker.loc[i, f'{factor_name}_start'] = overlap_row['start']
            overlap_tracker.loc[i, f'{factor_name}_end'] = overlap_row['end']
            overlap_tracker.loc[i, f'{factor_name}_center'] = overlap_row['center']
        else:
            print(f'Warning: for ARS_id {i} found more than one overlappin interval for {factor_name}. Saving only first (edit this in the future)')
            overlap_row = curr_overlap.iloc[0]

            overlap_tracker.loc[i, f'{factor_name}_start'] = overlap_row['start']
            overlap_tracker.loc[i, f'{factor_name}_end'] = overlap_row['end']
            overlap_tracker.loc[i, f'{factor_name}_center'] = overlap_row['center']
    
    return pd.concat([origin_df, overlap_tracker], axis=1)
