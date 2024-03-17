
#%%
# rearragne the predictions in the basecalled fragments in the correct order back to one npz file
# call this file basecalled_train_data_rearranged.npz

import os
import numpy as np

basecalling_dir = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_basecalled'

def recombine_basecalling(input_dir, output_file_path):
    # List all .npz files and sort them by the fragment number
    npz_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.npz')],
                       key=lambda x: int(x.split('_')[-1].split('.')[0]))

    # Initialize a list to hold the recombined basecalling data for each sample later
    recombined_data = []

    # Loop over each sample index
    for sample_idx in range(32):  # Assuming there are 32 samples as mentioned
        # Initialize a list to hold fragments for the current sample later 
        sample_fragments = []

        # Loop through each sorted npz file and extract the corresponding sample's basecalling
        for file_name in npz_files:
            file_path = os.path.join(input_dir, file_name)
            data = np.load(file_path)
            sample_fragments.append(data['basecalling'][sample_idx])

        # Concatenate fragments for the current sample and add to the recombined list
        sample_long_read = np.concatenate(sample_fragments, axis=0)
        recombined_data.append(sample_long_read)

    # Convert the list of recombined data to a numpy array with the desired shape
    recombined_array = np.array(recombined_data)

    # Save the recombined array to an npz file
    np.savez_compressed(output_file_path, recombined_basecalling=recombined_array)



