#%%
import numpy as np
import matplotlib.pyplot as plt
import os

""" Main script to run LONG READ pipeline """

# working directory
base_path = os.getcwd()

# Steps of the pipeline

#%%
# 1. Generate LONG READ data (60_000 time points) - save in folder LONG_READ_training_data in wd

from LONG_READ_1_RNA_data_generation import Generate_sim_signal, kmer_info, path_save_long_read

seq_len = 3000 # was 35 for 400 - 800 time points (benchamrking), increase
time_points = 60_000  
N_batch = 1 # Number of npz files that are generated
batch_size = 32

Generate_sim_signal(N_batch, batch_size , path_save_long_read , kmer_info, seq_len, time_points)

#%%
# 2. Split LONG READ into fragments à 1200 time points - save in folder LONG_READ_training_data_splitted
from LONG_READ_2_Split import split_and_save, npz_files, path_to_long_read, path_save

fragment_length = 1200

for npz_file in npz_files:
    full_path_to_file = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path_to_file, fragment_length)


#%%
# 3. Basecalling on fragments - save results in folder
from LONG_READ_3_Basecalling import predict_and_save_basecalling, input_dir, output_dir, model_path

    # Ensure the output directory exists

predict_and_save_basecalling(input_dir, output_dir, model_path)


#%%
# 4. Rearrange results of fragments back to LONG READ
from LONG_READ_4_Rearrange import rearrange_fragments_to_long_read, input_dir, path_save_rearranged, output_filename

base_path = os.getcwd()

rearrange_fragments_to_long_read(input_dir, path_save_rearranged, output_filename)

#%%

# 4. Rearrange results of fragments back to LONG READ
rearranged_npz = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_rearranged\basecalled_long_read_rearranged.npz"
data = np.load(rearranged_npz)
print(list(data.keys()))

long_read = data['long_read']
print(long_read)

# long read is saved with 
"""
['long_read']
[[[1.50911538e-02 5.19840181e-01 5.17301844e-04 6.47091586e-03
   4.58080381e-01]
  [6.93952339e-03 4.34740514e-01 2.79572268e-04 5.56761539e-03
   5.52472770e-01]
  [3.66408913e-03 3.83335382e-01 2.74935010e-04 6.44197734e-03
   6.06283545e-01]
  ...
  [1.68079495e-01 3.43627006e-01 2.08609685e-01 2.32912704e-01
   4.67710607e-02]
  [1.62946343e-01 3.31625313e-01 2.07538679e-01 2.28884533e-01
   6.90051466e-02]
  [1.52573749e-01 3.03563654e-01 2.01893032e-01 2.24403754e-01
   1.17565751e-01]]

 [[2.81222351e-02 5.32392692e-03 5.82746685e-01 3.03320307e-03
   3.80774051e-01]
  [1.60369612e-02 3.69668752e-03 5.52399814e-01 3.24702542e-03
   4.24619466e-01]
  [1.17807081e-02 2.47922959e-03 5.66019773e-01 4.02401946e-03
   4.15696263e-01]
  ...
"""







#%%
import numpy as np
base_path = os.getcwd()
def load_npz_file(file_path):
    """
    Load an NPZ file and print its keys.
    """
    with np.load(file_path) as data:
        print("Keys in the NPZ file:", data.keys())
        return {key: data[key] for key in data.keys()}

output_dir = os.path.join(base_path, 'LONG_READ_training_data_rearranged')
load_npz_file(output_dir)

data = load_npz_file(output_dir)
basecalling = data['basecalling']
original_data = data['original_data']

print(basecalling)
print(original_data)

#%%
from LONG_READ_5_Original_vs_Basecalling import compare_sequences, load_sequence, plot_mismatches
base_path = os.getcwd()

# Load sequences
original_seq_path = os.path.join(base_path, 'LONG_READ_training_data')
basecalled_seq_path = os.path.join(base_path, 'LONG_READ_data_rearranged')

# Assuming you're comparing the first .npz file in each directory
original_seq_file = next((f for f in os.listdir(original_seq_path) if f.endswith('.npz')), None)
basecalled_seq_file = next((f for f in os.listdir(basecalled_seq_path) if f.endswith('.npz')), None)

if original_seq_file and basecalled_seq_file:
    original_seq = load_sequence(os.path.join(original_seq_path, original_seq_file), 'rand_seq')
    basecalled_seq = load_sequence(os.path.join(basecalled_seq_path, basecalled_seq_file), 'basecalling')

    mismatches = compare_sequences(original_seq, basecalled_seq)
    plot_mismatches(original_seq, basecalled_seq, mismatches)
else:
    print("Required .npz files not found in the directories.")
# %%
