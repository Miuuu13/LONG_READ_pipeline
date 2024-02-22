#%%

#Steps to take:

#%%
# 1. Generate LONG READ data (6000 time points) - save in folder LONG_READ_training_data

from LONG_READ_1_RNA_data_generation import Generate_sim_signal, N_batch, batch_size , path_save, kmer_info, seq_len, time_points

Generate_sim_signal(N_batch, batch_size , path_save, kmer_info, seq_len, time_points)

#%%
# 2. Split LONG READ into fragments à 800 time points - save in folder LONG_READ_training_data_splitted
from LONG_READ_2_Split import split_and_save, npz_files, path_to_long_read
import os

for npz_file in npz_files:
    full_path = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path, fragment_length=800, overlap = 0)
    print(f"Processed and split: {npz_file}")

#%%
# 3. Basecalling on fragments - save results in folder
from LONG_READ_3_Basecalling import predict_and_save_basecalling, input_dir, output_dir, model_path

    # Ensure the output directory exists


predict_and_save_basecalling(input_dir, output_dir, model_path)

#%%
# 4. Rearrange results of fragments back to LONG READ
from LONG_READ_4_Rearrange import rearrange_fragments
# Number of batches to process
N_batch = 3  # Assuming there are 3 long reads as per your setup
overlap = 50  # Overlap between consecutive fragments

rearrange_fragments(N_batch, overlap)
#%%
# 5. Compare original LONG READ rand_seq with basecalling prediction
