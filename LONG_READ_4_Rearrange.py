
#%%
import numpy as np
import os

def rearrange_fragments_to_long_read(input_dir, path_save_rearranged, output_filename):
    """
    Rearrange the split fragments back to a long read and save the concatenated
    basecalled sequence along with the original data from the first fragment.
    
    Args:
    - input_dir: Directory containing the split fragment files for basecalling.
    - path_save_rearranged: Directory where the rearranged long read will be saved.
    - output_filename: Filename for the saved rearranged long read.
    """
    # Check if the output directory exists, if not, create it
    if not os.path.exists(path_save_rearranged):
        os.makedirs(path_save_rearranged)

    # Initialize lists to hold the basecalled data and extract original data
    basecalling_data = []
    original_data = None  # Will hold the original data from the first fragment

    # Get all basecalled fragment filenames and sort them
    fragment_files = [f for f in os.listdir(input_dir) if f.startswith('basecalled_train_data_0_fragment_') and f.endswith('.npz')]
    fragment_files_sorted = sorted(fragment_files, key=lambda x: int(x.split('_')[-1].split('.')[0]))

    # Load each fragment to concatenate basecalling and extract original data from the first fragment
    for i, file_name in enumerate(fragment_files_sorted):
        file_path = os.path.join(input_dir, file_name)
        with np.load(file_path) as data:
            # Concatenate basecalling predictions
            basecalling = data['basecalling']  # Extract basecalling predictions
            basecalling_data.append(basecalling)
            
            # Extract original data from the first fragment
            if i == 0:  # If it's the first fragment, extract original data
                original_data = data['original_data']

    if not basecalling_data:
        raise ValueError("No fragments were loaded. Check the input directory and filenames.")

    # Concatenate all basecalling predictions to form the long read basecalling sequence
    concatenated_basecalling = np.concatenate(basecalling_data, axis=0)

    # Save the rearranged long read with concatenated basecalling and original data from the first fragment
    output_file_path = os.path.join(path_save_rearranged, output_filename)
    np.savez(output_file_path, basecalling=concatenated_basecalling, original_data=original_data)

# Example usage
base_path = os.getcwd()
input_dir = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
path_save_rearranged = os.path.join(base_path, 'LONG_READ_training_data_rearranged')
output_filename = 'basecalled_long_read_rearranged.npz'

rearrange_fragments_to_long_read(input_dir, path_save_rearranged, output_filename)
























#%%
import numpy as np
import os

def rearrange_fragments_to_long_read(input_dir, path_save_rearranged, output_filename):
    """
    Rearrange the split fragments back to a long read.
    
    Args:
    - original_file_prefix: Prefix of the original file names before splitting.
    - num_fragments: Number of fragments each file was split into.
    - overlap: Number of overlapping points between consecutive fragments.
    
    Returns:
    - reconstructed_signal: The reconstructed signal array.
    - reconstructed_map: The reconstructed map (one-hot encoded) array.
    - original_rand_seq: The original random sequence from the first fragment.
    """
    # Check if the output directory exists, if not, create it
    if not os.path.exists(path_save_rearranged):
        os.makedirs(path_save_rearranged)

    # Initialize an empty list to hold the data
    fragments_data = []

    # Get all basecalled fragment filenames and sort them
    fragment_files = [f for f in os.listdir(path_basecalled) if f.startswith('basecalled_train_data_0_fragment_') and f.endswith('.npz')]
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

    output_file_path = os.path.join(path_save_rearranged, output_filename)
    np.savez(output_file_path, long_read=long_read)

# Directories
base_path = os.getcwd()
path_basecalled = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
path_save_rearranged = os.path.join(base_path, 'LONG_READ_training_data_rearranged')  




#%%
