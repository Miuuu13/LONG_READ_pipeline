#%%
import numpy as np
import os

""" 7. Compare """
path_original = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
path_basecalled = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"

""" 
Check the content of the two npz files after loading 
print the keys, data, and data type
"""
def print_npz_keys(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: \n {list(data.keys())}")
            for key in data.keys():
                print(f"Content of {key}:{data[key]}")
                print(f"Shape of {key}:{data[key].shape}")
                print(f"Data type of {key}: {type(data[key])} ")
                
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print_npz_keys(path_original)
print_npz_keys(path_basecalled)




# %%

""" Load the data """
#original data (1 sequences)
original_data = np.load(path_original)
original_seq = original_data['rand_seq']
print(original_seq)

#basecalled data (32 sequences)
basecalled_data = np.load(path_basecalled)
basecalled_seq = basecalled_data['seq1']
print(basecalled_seq)

""" NEXT: Decode the original seq from numbers to letters """
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'} 
original_seq_str = ''.join(index_to_base[base] for base in original_seq)

print(original_seq_str)
# %%

""" Use remove consecutive duplicates """

def remove_consecutive_duplicates(s):
    result = [s[0]]  # Start with the first char
    for char in s[1:]:  # Iterate over the string, starting at the second char
        if char != result[-1]:  # If the current char is not the same as the last in result
            result.append(char)  # Append the current char to result
    return ''.join(result)

print(remove_consecutive_duplicates(original_seq_str))

# %%


#%%

""" Remove here consecutive duplicates in rand_seq before decoding into letters """

import numpy as np

path_original = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
path_basecalled = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"

def remove_consecutive_duplicates(s):
    return ''.join(ch for i, ch in enumerate(s) if i == 0 or ch != s[i-1])

def print_npz_keys(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: \n {list(data.keys())}")
            for key in data.keys():
                print(f"Content of {key}:")
                print(data[key])
                print(f"Shape of {key}: {data[key].shape}")
                print(f"Data type of {key}: {type(data[key])} ")
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print_npz_keys(path_original)
print_npz_keys(path_basecalled)

""" Load the data """
#original data (1 sequences)
original_data = np.load(path_original)
original_seq = original_data['rand_seq']

#basecalled data (32 sequences)
basecalled_data = np.load(path_basecalled)
basecalled_seq = basecalled_data['seq1']

""" NEXT: Decode the original seq from numbers to letters """
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'} 
original_seq_str = ''.join(index_to_base[base] for base in original_seq)
original_seq_str = remove_consecutive_duplicates(original_seq_str)

print("Original sequence string without consecutive duplicates:", original_seq_str)
print("Length of original sequence string without consecutive duplicates:", len(original_seq_str))
#%%
""" save str seq original to npz"""
base_path = os.getcwd()
path_save_original_seq_str = base_path
np.savez(path_save_original_seq_str, original_seq_str=original_seq_str)

print(f"Original sequence string saved to: {path_save_original_seq_str}")
original_seq_str_loaded = np.load(original_seq_str)
original_seq_for_comparison = original_seq_str_loaded['']
