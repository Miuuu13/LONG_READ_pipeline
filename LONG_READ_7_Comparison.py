#%%
import numpy as np
import os

""" 7. Compare """

#%%

""" Remove here consecutive duplicates in rand_seq before decoding into letters """


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

""" save str seq original to npz"""
base_path = os.getcwd()
path_save_original_seq_str = base_path
np.savez(path_save_original_seq_str, original_seq_str=original_seq_str)

print(f"Original sequence string saved to: {path_save_original_seq_str}")
original_seq_str_loaded = np.load(original_seq_str)
original_seq_for_comparison = original_seq_str_loaded['original_seq_str']

#%%

""" Use pairwise for comparison """
import os
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def remove_consecutive_duplicates(s):
    return ''.join(ch for i, ch in enumerate(s) if i == 0 or ch != s[i-1])

"""
def print_original_vs_basecalled(original_seq_str, basecalled_data):
    print("Original sequence for comparison:")
    print(original_seq_str)
    print("Number of letters in original sequence:", len(original_seq_str))
    print("\nBasecalled sequences:")
    for key, seq in basecalled_data.items():
        print(f"Key: {key}")
        print("Sequence:", seq)
        print("Number of letters in sequence:", len(seq) if isinstance(seq, np.ndarray) else len(seq[()]))
        print()
"""
def calculate_alignment_stats(original_seq_str, basecalled_data):
    alignment_stats = {}
    for key, seq in basecalled_data.items():
        seq_str = str(seq)
        alignments = pairwise2.align.globalms(original_seq_str, seq_str, 2, -1, -0.5, -0.1)
        alignment = alignments[0]  # Considering the first alignment for simplicity
        alignment_str = format_alignment(*alignment)
        accuracy = (alignment[2] / len(original_seq_str)) * 100
        mismatches = alignment_str.count('| ')  # Mismatches represented as '|' in the alignment
        insertions = seq_str.count('-')  # Insertions represented as '-' in the basecalled sequence
        deletions = original_seq_str.count('-')  # Deletions represented as '-' in the original sequence
        alignment_stats[key] = {'Accuracy': accuracy, 'Mismatches': mismatches, 'Insertions': insertions, 'Deletions': deletions}
    return alignment_stats


# Paths
path_original = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
path_basecalled = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"

# Load data
original_data = np.load(path_original)
original_seq = original_data['rand_seq']
original_seq_str = ''.join(['ACGT'[base] for base in original_seq])
original_seq_str = remove_consecutive_duplicates(original_seq_str)

basecalled_data = np.load(path_basecalled)

# 1. Print original vs basecalled sequences
# print_original_vs_basecalled(original_seq_str, basecalled_data)

# 2. Calculate alignment stats
alignment_stats = calculate_alignment_stats(original_seq_str, basecalled_data)
for key, stats in alignment_stats.items():
    print(f"\nAlignment stats for {key}:")
    print(f"Accuracy: {stats['Accuracy']:.2f}%")
    print(f"Mismatches: {stats['Mismatches']}")
    print(f"Insertions: {stats['Insertions']}")
    print(f"Deletions: {stats['Deletions']}")











#%%


#%%





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#%%
