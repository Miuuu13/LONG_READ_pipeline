#%% 
import os
import numpy as np
#Code to split the long read
# long read data is saved in folder LONG_READ_training_data in working directory
# Plan: Split the long read into fragments à 800 time points with an overlap of 
base_path = os.getcwd()
#path_save = path_in + r"\Simulated_data_30JAN24_batch32\Simulated_data_n2_w800_b32"

path_to_long_read = os.path.join(base_path, 'LONG_READ_training_data')
path_save = os.path.join(base_path, 'LONG_READ_training_data_splitted')
#os.path.join(path_in, "Simulated_LONG_READ_data_FEB24", "Simulated_data_n2_w800_b32_split")

# Check if the save path exists, if not, create it
if not os.path.exists(path_save):
    os.makedirs(path_save)

# List all .npz files in the directory
npz_files = [f for f in os.listdir(path_to_long_read) if f.endswith('.npz')]
print(npz_files)

#%%

# use this function later to loop through all files in the directory that contains the long read data
def split_and_save(npz_file, segment_length, overlap):
    data = np.load(npz_file)
    signal_train = data['signal_train']
    map_onehot = data['map_onehot']
    rand_seq = data['rand_seq']
    
    num_segments = int((signal_train.shape[1] - segment_length) / (segment_length - overlap)) + 1
    
    base_name = npz_file[:-4]  # Remove .npz extension
    for i in range(num_segments):
        start = i * (segment_length - overlap)
        end = start + segment_length
        if end > signal_train.shape[1]:
            break  # Stop if the next segment would go beyond the array length
        
        segment_signal = signal_train[:, start:end]
        segment_map = map_onehot[:, start:end, :]
        #segment_rand_seq = rand_seq[:, start:end, :]
        
        # Save the segment
        segment_file_name = f"{base_name}_{chr(97 + i)}.npz"  # a, b, c, ...
        np.savez_compressed(segment_file_name, signal_train=segment_signal, map_onehot=segment_map, rand_seq=rand_seq)
        print(f"Saved: {segment_file_name}")


segment_length = 800
overlap = 50

# Process each file
for npz_file in npz_files:
    full_path = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path, segment_length, overlap)
    print(f"Processed and split: {npz_file}")







# %%
