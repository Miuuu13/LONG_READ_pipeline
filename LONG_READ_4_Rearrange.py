
#%%
import numpy as np
import os

def rearrange_fragments_to_long_read(input_dir, output_dir, output_filename):
    # Check if the output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize an empty list to hold the data
    fragments_data = []

    # Get all basecalled fragment filenames and sort them
    fragment_files = [f for f in os.listdir(input_dir) if f.startswith('basecalled_train_data_0_fragment_') and f.endswith('.npz')]
    fragment_files_sorted = sorted(fragment_files, key=lambda x: int(x.split('_')[-1].split('.')[0]))

    # Load each fragment and append its basecalling to the list
    for file_name in fragment_files_sorted:
        file_path = os.path.join(input_dir, file_name)
        with np.load(file_path) as data:
            basecalling = data['basecalling']  # Assuming 'basecalling' is the key for predictions
            fragments_data.append(basecalling)

    # Concatenate all fragments to form the long read
    long_read = np.concatenate(fragments_data, axis=0)

    # Save the rearranged long read
    output_file_path = os.path.join(output_dir, output_filename)
    np.savez(output_file_path, long_read=long_read)

# Directories
base_path = os.getcwd()
input_dir = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
output_dir = os.path.join(base_path, 'LONG_READ_training_data_rearranged')  # Updated output directory




#%%
import os
import numpy as np

# Setup the paths
base_path = os.getcwd()
split_data_path = os.path.join(base_path, 'LONG_READ_training_data_splitted')
rearranged_data_path = os.path.join(base_path, 'LONG_READ_rearranged')
original_data_path = os.path.join(base_path, 'LONG_READ_training_data')

# Check if the rearranged data path exists, if not, create it
if not os.path.exists(rearranged_data_path):
    os.makedirs(rearranged_data_path)

def rearrange_fragments(N_batch, overlap):
    for batch_id in range(N_batch):
        # Initialize containers for the reassembled signal and map
        reassembled_signal = np.array([])
        reassembled_map = np.array([])
        
        # Process each part
        for part in ['a', 'b', 'c', 'd', 'e', 'f', 'g']:
            file_name = f"basecalled_train_data_{batch_id}_{part}.npz"
            file_path = os.path.join(split_data_path, file_name)
            if os.path.exists(file_path):
                data = np.load(file_path)
                signal_fragment = data['signal_train']
                map_fragment = data['map_onehot']
                
                # Handle the overlap
                if reassembled_signal.size == 0:
                    reassembled_signal = signal_fragment
                    reassembled_map = map_fragment
                else:
                    # Adjust for overlap
                    overlap_adjustment = overlap // 2  # Adjust this based on how you want to handle the overlap
                    reassembled_signal = np.concatenate((reassembled_signal[:, :-overlap_adjustment], signal_fragment[:, overlap_adjustment:]), axis=1)
                    reassembled_map = np.concatenate((reassembled_map[:, :-overlap_adjustment, :], map_fragment[:, overlap_adjustment:, :]), axis=1)
        
        # Save the reassembled long read
        reassembled_file_name = f"rearranged_train_data_{batch_id}.npz"
        reassembled_file_path = os.path.join(rearranged_data_path, reassembled_file_name)
        np.savez_compressed(reassembled_file_path, signal_train=reassembled_signal, map_onehot=reassembled_map)
        print(f"Saved reassembled file: {reassembled_file_path}")








#%%
import numpy as np
import os

def rearrange_fragments(base_path, original_file_prefix, num_fragments=7, overlap=50):
    """
    Rearrange the split fragments back to a long read.
    
    Args:
    - base_path: Path to the directory containing the split .npz files.
    - original_file_prefix: Prefix of the original file names before splitting.
    - num_fragments: Number of fragments each file was split into.
    - overlap: Number of overlapping points between consecutive fragments.
    
    Returns:
    - reconstructed_signal: The reconstructed signal array.
    - reconstructed_map: The reconstructed map (one-hot encoded) array.
    - original_rand_seq: The original random sequence from the first fragment.
    """
    reconstructed_signal = []
    reconstructed_map = []
    
    for i in range(num_fragments):
        fragment_filename = f"{original_file_prefix}_{chr(97 + i)}.npz"  # a, b, c, ...
        fragment_path = os.path.join(base_path, fragment_filename)
        
        with np.load(fragment_path) as data:
            signal_train = data['signal_train']
            map_onehot = data['map_onehot']
            if i == 0:  # Only take rand_seq from the first fragment
                original_rand_seq = data['rand_seq']
            
            if i > 0:  # Adjust for overlap
                signal_train = signal_train[:, overlap:]
                map_onehot = map_onehot[:, overlap:, :]
            
            reconstructed_signal.append(signal_train)
            reconstructed_map.append(map_onehot)
    
    # Concatenate along the time axis
    reconstructed_signal = np.concatenate(reconstructed_signal, axis=1)
    reconstructed_map = np.concatenate(reconstructed_map, axis=1)
    
    return reconstructed_signal, reconstructed_map, original_rand_seq

# %%
