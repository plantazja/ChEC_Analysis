import os
import pandas as pd

def load_samples(*args):
    # Init directory to keeping samples df and corresponding names
    samples_dict = {}

    for s_path in args:
        if s_path is not None:
            # Check if samples are bed files and exists
            if not os.path.exists(s_path) and not s_path.endswith('.bed'):
                raise Exception(f'Provided file {s_path} does not exists or it is not bed file.')
            
            s_name = os.path.basename(s_path)
            samples_dict[s_name] = pd.read_csv(s_path, names=['chr', 'start', 'end', 'center'], sep='\t', header=None)
    return samples_dict
    