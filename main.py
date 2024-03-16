#%%
import numpy as np
import matplotlib.pyplot as plt
import os

""" Main script to run LONG READ pipeline (includes plots and checks that have to be removed later)"""
# working directory (wd)
base_path = os.getcwd()

# fixed parameters:
seq_len = 3000 # was 35 for 400 - 800 time points (benchamrking), increase for long read
time_points = 60_000  
N_batch = 1 # Number of npz files that are generated
batch_size = 32

#%%
""" 1. Generate LONG READ data (60_000 time points) - save in folder LONG_READ_training_data in wd """

from LONG_READ_1_Generate_RNA_Long_read import Generate_sim_signal, kmer_info, path_save_long_read

Generate_sim_signal(N_batch, batch_size , path_save_long_read , kmer_info, seq_len, time_points)

#%%
""" Plot the signal of the generated long read """
directory_path = path_save_long_read #folder 'LONG_READ_training_data'
file_name = 'train_data_0.npz'
file_path = os.path.join(directory_path, file_name)

data = np.load(file_path)
signal_train = data['signal_train'][0]  

def plot_signal_train(signal_train, file_name):
    plt.figure(figsize=(14, 5))  
    plt.plot(signal_train)
    plt.title(f'Signal Train for {file_name}')
    plt.xlabel('Time Points')
    plt.ylabel('Signal Value')
    plt.show()

plot_signal_train(signal_train, file_name)

#%%
def print_npz_contents(file_path):
    
    data = np.load(file_path)
    
    # list all keys
    print("Key in .npz file:")
    for key in data.keys():
        print(f"- {key}")
    
    # check content of each key
    for key in data.keys():
        print(f"\nContent of '{key}':")
        content = data[key]
        print(content)
        
        print(f"Data type: {content.dtype}, Shape: {content.shape}")

print_npz_contents(file_path)
#%%
""" 2. Split LONG READ into fragments à 1200 time points - save in folder LONG_READ_training_data_splitted """
from LONG_READ_2_Split import split_and_save, npz_files, path_to_long_read

fragment_length = 1200

for npz_file in npz_files:
    full_path_to_file = os.path.join(path_to_long_read, npz_file)
    split_and_save(full_path_to_file, fragment_length)

#%%
    
""" Check a splitted npz file """

def print_npz_contents(file_path):
    
    data = np.load(file_path)
    
    # list all keys
    print("Key in .npz file:")
    for key in data.keys():
        print(f"- {key}")
    
    # check content of each key
    for key in data.keys():
        print(f"\nContent of '{key}':")
        content = data[key]
        print(content)
        
        print(f"Data type: {content.dtype}, Shape: {content.shape}")
        

file_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_splitted\train_data_0_fragment_3.npz"
print_npz_contents(file_path)

#%%
""" Plot first split (fragment 1) signal of the long read"""  
 
from LONG_READ_plotting_functions import plot_signal_train 
directory_path = os.path.join(base_path, "LONG_READ_training_data_splitted")
file_name = 'train_data_0_fragment_1.npz'
file_path = os.path.join(directory_path, file_name)

data = np.load(file_path)
signal_train = data['signal_train'][0]  #i plotted several, to check

plot_signal_train(signal_train, file_name)

#%%
""" 3. Basecalling on fragments - save results in folder"""
from LONG_READ_3_Basecalling import predict_and_save_basecalling, input_dir, output_dir, model_path

predict_and_save_basecalling(input_dir, output_dir, model_path)

#%%

def print_npz_contents(file_path):
    
    data = np.load(file_path)
    
    # list all keys
    print("Key in .npz file:")
    for key in data.keys():
        print(f"- {key}")
    
    # check content of each key
    for key in data.keys():
        print(f"\nContent of '{key}':")
        content = data[key]
        print(content)
        
        print(f"Data type: {content.dtype}, Shape: {content.shape}")

file_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_basecalled\basecalled_train_data_0_fragment_1.npz"
print_npz_contents(file_path)

#%%

"""check basecalling - looks weird """


#TODO: Change to contour plot or heatmap (after integrating step 6)

path_basecalled_frament = r"LONG_READ_training_data_basecalled/basecalled_train_data_0_fragment_10.npz"
data = np.load(path_basecalled_frament)

predictions = data['basecalling']

# Extend the base dictionary to include also the spacer
base_dict = {0: "A", 1: "C", 2: "G", 3: "T", 4: "Spacer"}
# Define colors for plotting, adding an additional color for the spacer (grey or black?)
colors = ['blue', 'green', 'red', 'yellow', 'gray']  

# Plot the predictions for each base including the spacer
for base_index, base_label in base_dict.items():
    plt.plot(predictions[:, base_index], label=base_label, color=colors[base_index], linestyle='-', marker='')

plt.title('Basecalling Predictions for the First Fragment')
plt.xlabel('Sequence Position')
plt.ylabel('Prediction Score')

# Add the legend outside the loop to avoid duplicate legend entries
plt.legend(title="Base")
plt.show()


#%%

""" 4. Rearrange results of fragments back to LONG READ """
from LONG_READ_4_Rearrange import recombine_basecalling

# Use the function to recombine basecalling data
base_path = os.getcwd()
input_dir = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
output_file_path = os.path.join(base_path, 'recombined_long_read_basecalling.npz')

recombine_basecalling(input_dir, output_file_path)

#%%

"""Check recombination and the shape"""
path_recombined = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\recombined_long_read_basecalling.npz"

basecalling_data = np.load(path_recombined)['recombined_basecalling']

def print_npz_contents(file_path):
    
    data = np.load(file_path)
    
    # list all keys
    print("Key in .npz file:")
    for key in data.keys():
        print(f"- {key}")
    
    # check content of each key
    for key in data.keys():
        print(f"\nContent of '{key}':")
        content = data[key]
        print(content)
        
        print(f"Data type: {content.dtype}, Shape: {content.shape}")

print_npz_contents(path_recombined)

#%%

""" 5. Get basecalled sequence using argmax - printing the seq"""

import numpy as np

path_recombined = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\recombined_long_read_basecalling.npz"

# Load the .npz file and extract the 'recombined_basecalling' data
basecalling_data = np.load(path_recombined)['recombined_basecalling']

# Dictionary to map indices to nucleotide bases, including the spacer 'S'
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}

# Function to decode numerical array to sequence, including the spacer 'S'
def decode_sequence(numerical_array):
    return ''.join(index_to_base[idx] for idx in numerical_array)

# Decode sequences for all 32 samples
decoded_sequences = []
for sample in basecalling_data:
    # Use argmax to find the index of the highest value at each time point
    predicted_indices = np.argmax(sample, axis=1)
    # Decode the indices to a sequence
    sequence = decode_sequence(predicted_indices)
    decoded_sequences.append(sequence)

# Print the decoded sequences for all 32 samples
for i, seq in enumerate(decoded_sequences):
    print(f"Sample {i+1} Sequence: {seq}")


#%%
    
""" 5. Get basecalled sequence using argmax - saving the 32seqs into a npz file -
This was maily for checking how the sequences now look like the next cell is the one that saves them into a npz file 
before doing the comparison, there need to be removed: consecutive nases and the spacer!"""
import numpy as np

# Full path to your .npz file
path_recombined = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\recombined_long_read_basecalling.npz"

basecalling_data = np.load(path_recombined)['recombined_basecalling']

# map indices to nucleotide bases, including the spacer 'S' using a dict.
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}

# Function to decode numerical array to sequence
def decode_sequence(numerical_array):
    return ''.join(index_to_base[idx] for idx in numerical_array)

# Decode sequences for all 32 samples
decoded_sequences = {}
for sample_idx, sample in enumerate(basecalling_data):
    # Use argmax to find the index of the highest value at each time point
    predicted_indices = np.argmax(sample, axis=1)
    # Decode the indices to a sequence
    sequence = decode_sequence(predicted_indices)
    # Store each sequence in the dictionary with a unique key name 
    decoded_sequences[f'seq{sample_idx+1}'] = sequence

# Save seqs into npz
""" check the keys of the seqs """
output_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
np.savez_compressed(output_path, **decoded_sequences)



#%%
npz_file_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
npz_file = np.load(npz_file_path)

# List all entries/keys in the npz file
keys = npz_file.files
print(keys)
first_key_shape = npz_file['seq1'].shape



#%%
import numpy as np 

""" Remove consecutive duplicated in the sequences"""

import numpy as np

def remove_consecutive_duplicates(s):
    result = [s[0]]  # Start with the first char
    for char in s[1:]:  # Iterate over the string, starting at the second char
        if char != result[-1]:  # If the current char is not the same as the last in result
            result.append(char)  # Append the current char to result
    return ''.join(result)

path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

shortened_sequences = {}


for key in npz_file.files:
    seq = npz_file[key]
    if isinstance(seq, np.ndarray) and seq.shape == ():  # Convert scalar arrays to Python scalars
        seq = seq.item()
    # remove the unwanted consecutive duplicates
    shortened_seq = remove_consecutive_duplicates(seq)
    shortened_sequences[key] = shortened_seq

output_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences.npz"
np.savez_compressed(output_path, **shortened_sequences)


#%%

""" Remove the spacer"""
import numpy as np

def remove_consecutive_duplicates_and_s(s):
    result = []  # Start with empty list
    previous_char = None  # Keep track of the previous char
    
    for char in s:
        # Check if the current char is different from the previous char and is not 'S'
        if char != previous_char and char != 'S':
            result.append(char)  # Append current char to the result
            previous_char = char  # Update the previous char
            
    return ''.join(result)

path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

shortened_sequences_without_S = {}


for key in npz_file.files:
    seq = npz_file[key]
    if isinstance(seq, np.ndarray) and seq.shape == ():  # Convert scalar arrays to Python scalars
        seq = seq.item()
    
    shortened_seq_without_S = remove_consecutive_duplicates_and_s(seq)
    shortened_sequences_without_S[key] = shortened_seq_without_S

output_path_without_S = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"
np.savez_compressed(output_path_without_S, **shortened_sequences_without_S)


#%%

#----------------- new 16MRC24 ----------------------------

""" 7. Compare """
path_original = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
path_basecalled = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"


def print_npz_keys(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: {list(data.keys())}")
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print_npz_keys(path_original)
print_npz_keys(path_basecalled)

#%%
def print_npz_keys_and_shapes(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: {list(data.keys())}")
            for key in data.keys():
                print(f"Shape of {key}: {data[key].shape}")
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print(f"Shape of rand_seq in {path_original}:")
with np.load(path_original) as data_original:
    print(data_original["rand_seq"].shape)

print_npz_keys_and_shapes(path_basecalled)

#%%
def print_npz_keys_and_content(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: {list(data.keys())}")
            for key in data.keys():
                print(f"Content of {key}:")
                print(data[key])
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print(f"Content of rand_seq in {path_original}:")
with np.load(path_original) as data_original:
    print(data_original["rand_seq"])

print_npz_keys_and_content(path_basecalled)

#%%

def print_npz_keys_and_content_length(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: {list(data.keys())}")
            for key in data.keys():
                print(f"Content of {key}:")
                print(data[key])
                print(f"Length of {key}: {len(data[key])}")
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print(f"Content of rand_seq in {path_original}:")
with np.load(path_original) as data_original:
    print(data_original["rand_seq"])

print_npz_keys_and_content_length(path_basecalled)



#%%
def print_npz_keys_and_content_length(file_path):
    try:
        with np.load(file_path) as data:
            print(f"Keys in {file_path}: {list(data.keys())}")
            for key in data.keys():
                print(f"Content of {key}:")
                print(data[key])
                if isinstance(data[key], np.ndarray):
                    print(f"Length of {key}: {data[key].shape}")
                else:
                    print(f"Length of {key}: Not applicable (not an ndarray)")
    except IOError:
        print(f"Error: Could not load data from {file_path}")

print(f"Content of rand_seq in {path_original}:")
with np.load(path_original) as data_original:
    print(data_original["rand_seq"])

print_npz_keys_and_content_length(path_basecalled)













#%%




import numpy as np
""" 6. Compare original vs basecalled sequence """
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']
print(original_seq)

path_rearranged = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
npz_file = np.load(path_rearranged)
keys = npz_file.files

print(keys)


#%%

import numpy as np

train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['signal_train']

path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

for key in npz_file.files:
    seq = npz_file[key]
    print("Original Seq:", original_seq)
    print(f"{key}:", seq)

#%%
import numpy as np

train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']

index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'} 
original_seq_str = ''.join(index_to_base[base] for base in original_seq)

path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)


seq1 = npz_file['seq1']
if isinstance(seq1, np.ndarray) and seq1.shape == ():  
    seq1 = seq1.item()
print("Length of seq1:", len(seq1))

print("Shape of the original sequence:", original_seq.shape)

#%%
import numpy as np

train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']

index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
original_seq_str = ''.join(index_to_base.get(base, '') for base in original_seq)  # Use .get() to handle missing keys

path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"
npz_file = np.load(path_rearranged)

seq1 = npz_file['seq1']
if isinstance(seq1, np.ndarray) and seq1.shape == (): 
    seq1 = seq1.item()
print("Length of seq1 (after processing):", len(seq1))

print("Length of the original sequence string:", len(original_seq_str))


#%%
""" 6. Compare original vs basecalled sequence using levenstein """
import numpy as np
from Levenshtein import distance as levenshtein_distance

def decode_rand_seq(rand_seq, index_to_base):
    return ''.join(index_to_base[base] for base in rand_seq)

def compare_sequences(original_sequence, decoded_sequences_file):
    comparison_results = {}
    for key in decoded_sequences_file.files:
        decoded_sequence = decoded_sequences_file[key]
        if isinstance(decoded_sequence, np.ndarray) and decoded_sequence.shape == ():
            decoded_sequence = decoded_sequence.item()  
        
        distance = levenshtein_distance(original_sequence, decoded_sequence)
        comparison_results[key] = distance
    return comparison_results

decoded_sequences_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
decoded_sequences_file = np.load(decoded_sequences_path)

train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
rand_seq = train_data_file['rand_seq']
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}
original_sequence = decode_rand_seq(rand_seq, index_to_base)

comparison_results = compare_sequences(original_sequence, decoded_sequences_file)

for seq_key, distance in comparison_results.items():
    print(f"Edit distance for {seq_key}: {distance}")

#%%%

import numpy as np
from Levenshtein import distance as levenshtein_distance

decoded_sequences_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\recombined_long_read_basecalling.npz"
decoded_sequences_file = np.load(decoded_sequences_path)

train_data_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\LONG_READ_training_data\\train_data_0.npz"
train_data_file = np.load(train_data_path)
rand_seq = train_data_file['signal_train']

index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}
original_sequence = ''.join(index_to_base[base] for base in rand_seq)
print(original_sequence)
def compare_sequences(decoded_sequences_file, original_sequence):
    results = {}
    for key in decoded_sequences_file.files:
        decoded_sequence = decoded_sequences_file[key]
        if isinstance(decoded_sequence, np.ndarray) and decoded_sequence.shape == ():
            decoded_sequence = decoded_sequence.item()  
            
        distance = levenshtein_distance(decoded_sequence, original_sequence)
        results[key] = distance
    return results

comparison_results = compare_sequences(decoded_sequences_file, original_sequence)

for seq_key, distance in comparison_results.items():
    print(f"Edit distance for {seq_key}: {distance}")

# %%

""" 6. Compare original vs basecalled sequence using pairwise """

import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

decoded_sequences_path = r"C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
decoded_sequences_file = np.load(decoded_sequences_path)

train_data_path = r"C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\LONG_READ_training_data\\train_data_0.npz"
train_data_file = np.load(train_data_path)
rand_seq = train_data_file['signal_train']
print(rand_seq)

index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}  
original_sequence = ''.join(index_to_base[base] for base in rand_seq)

def compare_sequences_pairwise(decoded_sequences_file, original_sequence):
    alignments = {}
    for key in decoded_sequences_file.files:
        decoded_sequence = decoded_sequences_file[key]
        if isinstance(decoded_sequence, np.ndarray) and decoded_sequence.shape == ():
            decoded_sequence = decoded_sequence.item()  

        alignment = pairwise2.align.globalxx(original_sequence, decoded_sequence, one_alignment_only=True)[0]
        alignments[key] = alignment
    return alignments

alignments = compare_sequences_pairwise(decoded_sequences_file, original_sequence)

for seq_key, alignment in alignments.items():
    print(f"Alignment for {seq_key}:")
    print(format_alignment(*alignment))
    print("\n")



#%%

""" How to force installation inside an environment """
#!pip install biopython


# %%
