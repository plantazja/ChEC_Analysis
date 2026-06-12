import os
import pandas as pd
import numpy as np
import itertools


def save_overlaps(overlaps_df, output_dir):
    sample_names = list(overlaps_df.columns[4:]) 
    sample_names = [n.replace('.bed', '') for n in sample_names]
    sample_indices = [i for i in range(len(sample_names))]
    
    # Create output subdirectory for given samples
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    str_names = '_'.join(sample_names)
    output_dir = os.path.join(output_dir, str_names)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Saving the files
    decoding = {} # Init dict: {code : full_name}
    code_dfs = {} # Init dict: {full_name : df}
    
    # 1. Find all possible combinations
    combs = []
    for i in range(1, len(sample_indices)+1):
        new_combs  = [list(c) for c in itertools.combinations(sample_indices, i)] # generate all combinations of length i
        combs.append(new_combs)
    
    # Sort each list of combinations and remove duplicates
    all_combs = []
    for comb_list in combs:
        for comb in comb_list:
            all_combs.append(tuple(sorted(comb, reverse=True)))

    combs_sorted = set(all_combs)

    # 2. Create decoding dict to go fast from code to full name ex.'in_sample1_ex_sample2'
    for comb in combs_sorted:
        # Init code for new combination
        code = [0] * len(sample_names)
        full_name_in = ['in']
        full_name_ex = ['ex']
        included_samples = [sample_names[i] for i in comb]
        excluded_samples = [sample for sample in sample_names if sample not in included_samples]
        
        # Update code and add samples to "include" list names
        for c in comb:
            code[c] = 1
        full_name_in.extend(included_samples)
        full_name_ex.extend(excluded_samples) # add samples to "exclude" that are not in combination
        
        full_name = full_name_in + full_name_ex
        full_name = '_'.join(full_name)

        # Convert code to tuple for use as key (lists are not hashable)
        decoding[tuple(code)] = full_name
        # Init dict for each full_name, in the end converting dict to pd.DataFrame
        code_dfs[full_name] = {'chr': [], 'start': [], 'end': [], 'center': []}
    
    # 3. Now iterate through overlap_df rows
    for _, row in overlaps_df.iterrows():
        row_code = tuple(row[4:].astype(int).tolist())
        if row_code in decoding:
            row_name = decoding[row_code]

            code_dfs[row_name]['chr'].append(row['chr'])
            code_dfs[row_name]['start'].append(row['start'])
            code_dfs[row_name]['end'].append(row['end'])
            code_dfs[row_name]['center'].append(row['center'])
    
    # 4. Save all created df with corresponding names
    #    also init dict to store peaks count for each combination
    dict_peaks = {'sample': [], 'peaks': []}
    for df_name, df_dict in code_dfs.items():
        df = pd.DataFrame(df_dict)
        file_name = df_name + '.bed'
        file_path = os.path.join(output_dir, file_name)
        df.to_csv(file_path, sep='\t', index=False, header=False)

        dict_peaks['sample'].append(df_name)
        dict_peaks['peaks'].append(df.shape[0])
    
    pd.DataFrame(dict_peaks).to_csv(os.path.join(output_dir, 'peaks_number.csv'), index=False)
    return output_dir