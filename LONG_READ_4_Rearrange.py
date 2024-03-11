
#%%
# rearragne the predictions in the basecalled fragments in the correct order back to one npz file
# call this file basecalled_train_data_rearranged.npz

import os
import numpy as np


basecalling_dir = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_basecalled'



def rearrange_basecalled_fragments(basecalling_dir):
    """
    Combines basecalling information from multiple .npz files in the correct order.

    Args:
        basecalling_dir (str): Path to the directory containing .npz files.

    Returns:
        None: Saves the rearranged basecalling data to 'rearranged_long_read.npz'.
    """
    # List all .npz files in the directory
    npz_files = [file for file in os.listdir(basecalling_dir) if file.endswith('.npz')]

    # Sort the files by fragment number
    sorted_npz_files = sorted(npz_files, key=lambda x: int(x.split('_')[-1].split('.')[0]))

    combined_basecalling = []
    rand_seq = None

    # Loop through all .npz files
    for npz_file in sorted_npz_files:
        file_path = os.path.join(basecalling_dir, npz_file)
        data = np.load(file_path)

        # Extract basecalling info
        basecalling_array = data['basecalling']
        print(basecalling_array)
        combined_basecalling.append(basecalling_array)
        

        # Extract rand_seq once (same for all fragments)
        if rand_seq is None:
            rand_seq = data['rand_seq']

    # Combine basecalling arrays in the correct order
    combined_basecalling = np.concatenate(combined_basecalling, axis=0)

    # Save the combined data to a new .npz file
    output_file = 'rearranged_long_read.npz'
    np.savez(output_file, basecalling=combined_basecalling, rand_seq=rand_seq)

    print(f"The rearranged basecalling information is now saved as '{output_file}'.")

# Example usage:
rearrange_basecalled_fragments(basecalling_dir)





#%%
import numpy as np
import matplotlib.pyplot as plt

data = np.load("rearranged_long_read.npz")

for key in data.keys():
    print(f"Key: {key}")
    print(f"Data type: {data[key].dtype}")
    print(f"Shape: {data[key].shape}\n")

#%%


