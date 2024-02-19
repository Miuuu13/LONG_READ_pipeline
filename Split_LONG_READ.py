#%%
import numpy as np
import os

base_path = os.getcwd()
#path_save = path_in + r"\Simulated_data_30JAN24_batch32\Simulated_data_n2_w800_b32"
path_in = os.path.join(base_path, 'LONG_READ_training_data')
# New directory for saving the fragments
split_dir = os.path.join(base_path, 'LONG_READ_training_data_splitted')

if not os.path.exists(split_dir):
    os.makedirs(split_dir)

def split_save_fragments(file_path, split_dir, window_size=800, overlap=750):
    # Load the npz file
    data = np.load(file_path)
    signal_train = data['signal_train']
    map_onehot = data['map_onehot']
    rand_seq = data['rand_seq'] if 'rand_seq' in data else None

    # Calculate number of fragments
    num_fragments = (signal_train.shape[1] - window_size) // (window_size - overlap) + 1

    # Extract the base file name without extension
    base_file_name = os.path.basename(file_path).replace('.npz', '')

    for i in range(num_fragments):
        start = i * (window_size - overlap)
        end = start + window_size
        signal_fragment = signal_train[:, start:end]
        map_fragment = map_onehot[:, start:end] if map_onehot is not None else None
        rand_fragment = rand_seq[start:end] if rand_seq is not None else None
        
        # Construct fragment file name
        fragment_file_name = f"{base_file_name}_{chr(97 + i)}.npz"
        
        # Save the fragment
        np.savez_compressed(os.path.join(split_dir, fragment_file_name),
                            signal_train=signal_fragment,
                            map_onehot=map_fragment,
                            rand_seq=rand_fragment)

# Process each of the npz files
for i in range(3):  # Assuming you have files named train_data_0.npz, train_data_1.npz, train_data_2.npz
    file_name = f"train_data_{i}.npz"
    file_path = os.path.join(base_path, file_name)
    split_save_fragments(file_path, split_dir)



