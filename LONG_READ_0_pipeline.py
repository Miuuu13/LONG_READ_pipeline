#%%
import numpy as np
import matplotlib.pyplot as plt
import os
# working directory
base_path = os.getcwd()

#Steps to take:

#%%
# 1. Generate LONG READ data (6000 time points) - save in folder LONG_READ_training_data

from LONG_READ_1_RNA_data_generation import Generate_sim_signal, N_batch, batch_size , path_save, kmer_info, seq_len, time_points

Generate_sim_signal(N_batch, batch_size , path_save, kmer_info, seq_len, time_points)

#%%
"""OPTIONAL: Plot
from Plot_npz_files import plot_all_data_in_directory

# directory containing my .npz files
directory_path = os.path.join(base_path, 'LONG_READ_training_data')

plot_all_data_in_directory(directory_path)
"""

#%%
# 2. Split LONG READ into fragments à 800 time points - save in folder LONG_READ_training_data_splitted
from LONG_READ_2_Split import split_and_save, npz_files, path_to_long_read
import os

for npz_file in npz_files:
    full_path = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path, fragment_length=1200, overlap = 0)
    print(f"Processed and split: {npz_file}")

#%%
# 3. Basecalling on fragments - save results in folder
from LONG_READ_3_Basecalling import predict_and_save_basecalling, input_dir, output_dir, model_path

    # Ensure the output directory exists

predict_and_save_basecalling(input_dir, output_dir, model_path)

#%%
from LONG_READ_4_Rearrange import  rearrange_fragments_to_long_read
# 4. Rearrange results of fragments back to LONG READ
# Call the function with the directory and the desired output filename
rearrange_fragments_to_long_read(input_dir, output_dir, 'basecalled_long_read_rearranged.npz')

#%%
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_data_rearranged\basecalled_long_read_rearranged.npz"



def plot_data(signal_train, map_onehot, rand_seq_numeric, file_name):
    # Convert numeric sequence to string using a predefined dictionary
    base_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    rand_seq = ''.join([base_dict[num] for num in rand_seq_numeric])
    
    # Start plotting
    fig, axs = plt.subplots(3, 1, figsize=(14, 20), constrained_layout=True)

    # Plot Signal Train
    axs[0].plot(signal_train)
    axs[0].set_title('Signal Train')

    # Plot Filled Contour for map_onehot
    if map_onehot.any():
        axs[1].contourf(map_onehot.T, cmap='viridis', levels=np.linspace(0, 1, num=50))
        axs[1].set_title('Filled Contour of One-hot Encoded Map')
        axs[1].set_xlabel('Time Point')
        axs[1].set_ylabel('Nucleotide Position')

    # Plot Heatmap for map_onehot using Seaborn
    sns.heatmap(map_onehot.T, cmap="YlGnBu", cbar_kws={'label': 'Feature Activation'}, ax=axs[2])
    axs[2].set_title('Map Onehot Features Over Time')
    axs[2].set_xlabel('Time Points')
    axs[2].set_ylabel('Features')
    # Adjusting the y-ticks to show all feature labels
    axs[2].set_yticks(range(map_onehot.shape[1]))
    axs[2].set_yticklabels([f'Feature {i}' for i in range(map_onehot.shape[1])])

    plt.suptitle(file_name)
    plt.show()

    # Print the sequence as a string
    print(f"Random Sequence for {file_name}: {rand_seq}")
    
#%%
# 5. Compare original LONG READ rand_seq with basecalling prediction
# %%
