import os

def save_signals(signals, in_file_path, out_dir_path):
    # Get file name
    filename = os.path.basename(in_file_path)
    sample = '_'.join(filename.split('_')[:6]) + '.csv'

    # Return output path for new .csv file
    out_path = os.path.join(out_dir_path, sample)
    os.makedirs(out_dir_path, exist_ok=True)

    # Save csv
    signals.to_csv(out_path, sep='\t', header=False, index=False)