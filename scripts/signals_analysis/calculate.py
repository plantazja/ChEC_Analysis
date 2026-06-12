import pandas as pd
import numpy as np

def get_signals(peaks_coord, signal_coord):
    CHROM = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX',
             'chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']
    
    all_peak_profiles = []
    
    for chr in CHROM:
        # Filter data for current chromosome
        chr_signals = signal_coord[signal_coord['chr'] == chr].copy()
        chr_peaks = peaks_coord[peaks_coord['chr'] == chr].copy()
        
        if chr_peaks.empty:
            continue
        
        # Pre-calculate intervals for faster lookup
        signal_starts = chr_signals['start'].values
        signal_ends = chr_signals['end'].values
        signal_values = chr_signals['signal'].values
        
        for _, row in chr_peaks.iterrows():
            start = row['start']
            end = row['end']
            peak_profile = []
            
            # For each position in the peak
            for pos in range(start, end + 1):
                # Search for intervals containing position
                mask = (signal_starts <= pos) & (signal_ends >= pos)
                
                if np.any(mask):
                    # Get first matching signal value
                    first_match_idx = np.where(mask)[0][0]
                    peak_profile.append(signal_values[first_match_idx])
                else:
                    # No interval contains this position - naive approach
                    peak_profile.append(0)
                    
            all_peak_profiles.append(peak_profile)
    
    signal_df = pd.DataFrame(all_peak_profiles).T
    if signal_df.shape[0] != 151:
        print('Peaks width is not 151 bp')
        signal_df = pd.DataFrame()
    return signal_df

def log2fc(IAA_df, DMSO_df, pseudo_count=0.0001):
    # Add pseudo count and log(0)
    IAA = IAA_df.sum(axis=0) + pseudo_count
    DMSO = DMSO_df.sum(axis=0) + pseudo_count 

    # Calculate fold change and log2 transform
    fold_change = IAA / DMSO 
    return np.log2(fold_change)