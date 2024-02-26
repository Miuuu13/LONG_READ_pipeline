#%%
import os
import numpy as np

base_path = os.getcwd()

path_to_long_read = os.path.join(base_path, r'LONG_READ_training_data')
path_save = os.path.join(base_path, r'LONG_READ_training_data_splitted')

# Check if the save path exists, if not, create it
if not os.path.exists(path_save):
    os.makedirs(path_save)

# List all .npz files in the directory
npz_files = [f for f in os.listdir(path_to_long_read) if f.endswith('.npz')]
print(npz_files)

def split_and_save(npz_file, fragment_length=1200, overlap=0):  # Default overlap to 0
    data = np.load(npz_file)
    signal_train = data['signal_train']
    map_onehot = data['map_onehot']
    rand_seq = data['rand_seq']
    
    # Update the calculation of num_fragments for no overlap
    num_fragments = int(np.ceil(signal_train.shape[1] / fragment_length))
    
    base_name = os.path.basename(npz_file)[:-4]  # Remove .npz extension and path
    for i in range(num_fragments):
        start = i * fragment_length
        end = start + fragment_length
        
        # Check if the end of the fragment exceeds the signal length
        if end > signal_train.shape[1]:
            print("The last fragment does not have 1200 time points. Stopping.")
            break  # Stop the loop if there are not enough time points left for a full fragment
        
        fragment_signal = signal_train[:, start:end]
        fragment_map = map_onehot[:, start:end, :]
        
        # Save the fragment with numerical suffix
        fragment_file_name = f"{base_name}_fragment_{i+1}.npz"  # Using numbers for fragment naming
        fragment_file_path = os.path.join(path_save, fragment_file_name)  # Path to save the fragment
        np.savez_compressed(fragment_file_path, signal_train=fragment_signal, map_onehot=fragment_map, rand_seq=rand_seq)
        print(f"Saved: {fragment_file_path}")

# Example usage, assuming you have the variables fragment_length and npz_file set appropriately.
# You would loop through `npz_files` and call `split_and_save` for each, like this:

fragment_length = 1200  # Example fragment length, adjust as needed

for npz_file in npz_files:
    full_path_to_file = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path_to_file, fragment_length)


