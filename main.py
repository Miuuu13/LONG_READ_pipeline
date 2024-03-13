#%%
import numpy as np
import matplotlib.pyplot as plt
import os

""" Main script to run LONG READ pipeline """

# working directory
base_path = os.getcwd()

# Steps of the pipeline

#%%
""" 1. Generate LONG READ data (60_000 time points) - save in folder LONG_READ_training_data in wd """

from LONG_READ_1_Generate_RNA_Long_read import Generate_sim_signal, kmer_info, path_save_long_read

seq_len = 3000 # was 35 for 400 - 800 time points (benchamrking), increase
time_points = 60_000  
N_batch = 1 # Number of npz files that are generated
batch_size = 32

Generate_sim_signal(N_batch, batch_size , path_save_long_read , kmer_info, seq_len, time_points)

#%%
""" Plot the signal of the generated long read """
#from LONG_READ_plotting_functions import plot_signal_train 
#import does not work...

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
# 2. Split LONG READ into fragments à 1200 time points - save in folder LONG_READ_training_data_splitted
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
        

# function call
#file_path = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_splitted\train_data_0_fragment_1.npz'
#print_npz_contents(file_path)
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
# 3. Basecalling on fragments - save results in folder
from LONG_READ_3_Basecalling import predict_and_save_basecalling, input_dir, output_dir, model_path

    # Ensure the output directory exists

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
        

# function call
#file_path = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_splitted\train_data_0_fragment_1.npz'
#print_npz_contents(file_path)
file_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_basecalled\basecalled_train_data_0_fragment_1.npz"
print_npz_contents(file_path)

#%%

# check basecalling
path_basecalled_frament = r"LONG_READ_training_data_basecalled/basecalled_train_data_0_fragment_10.npz"
data = np.load(path_basecalled_frament)

# 'basecalling' contains predictions with a shape (sequence_length, 5)
predictions = data['basecalling']

# Extend the base dictionary to include also the spacer
base_dict = {0: "A", 1: "C", 2: "G", 3: "T", 4: "Spacer"}
# Define colors for plotting, adding an additional color for the spacer (grey or black?)
colors = ['blue', 'green', 'red', 'yellow', 'gray']  # Assigning a unique color to each base and the spacer

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

# 4. Rearrange results of fragments back to LONG READ
from LONG_READ_4_Rearrange import recombine_basecalling

# Use the function to recombine basecalling data
base_path = os.getcwd()
input_dir = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
output_file_path = os.path.join(base_path, 'recombined_long_read_basecalling.npz')

recombine_basecalling(input_dir, output_file_path)


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

"""Check recombination and the shape"""
path_recombined = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\recombined_long_read_basecalling.npz"

# Load the .npz file, then access the 'recombined_basecalling' array
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

"""prediction plot ?"""
#%%

""" 5. Get basecalled sequence using argmax - printing the seq"""

import numpy as np

# Full path to your .npz file
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
    
""" 5. Get basecalled sequence using argmax - saving the 32seqs into a npz file """
import numpy as np

# Full path to your .npz file
path_recombined = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\recombined_long_read_basecalling.npz"

# Load the .npz file and extract the 'recombined_basecalling' data
basecalling_data = np.load(path_recombined)['recombined_basecalling']

# Dictionary to map indices to nucleotide bases, including the spacer 'S'
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}

# Function to decode numerical array to sequence, including the spacer 'S'
def decode_sequence(numerical_array):
    return ''.join(index_to_base[idx] for idx in numerical_array)

# Decode sequences for all 32 samples
decoded_sequences = {}
for sample_idx, sample in enumerate(basecalling_data):
    # Use argmax to find the index of the highest value at each time point
    predicted_indices = np.argmax(sample, axis=1)
    # Decode the indices to a sequence
    sequence = decode_sequence(predicted_indices)
    # Store each sequence in the dictionary with a unique key
    decoded_sequences[f'seq{sample_idx+1}'] = sequence

# Save the sequences into an npz file
output_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
np.savez_compressed(output_path, **decoded_sequences)



#%%
npz_file_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
npz_file = np.load(npz_file_path)

# List all entries/keys in the npz file
keys = npz_file.files
print(keys)
first_key_shape = npz_file['seq1'].shape
print(first_key_shape)

#empty shape
#%%
import numpy as np

# Load the .npz file
npz_file_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
npz_file = np.load(npz_file_path)

# Function to safely get the length of data under a key
def get_data_length(data):
    if isinstance(data, np.ndarray) and data.shape == ():  # If data is a scalar numpy array
        data = data.item()  # Convert to a Python scalar
    if isinstance(data, str):  # Check if the data is a string
        return len(data)
    else:
        return "Unknown format"  # Replace with appropriate handling

# Retrieve and print lengths for specific keys
first_key_length = get_data_length(npz_file['seq1'])
tenth_key_length = get_data_length(npz_file['seq10'])
last_key_name = sorted(npz_file.files)[-1]
last_key_length = get_data_length(npz_file[last_key_name])

print(f"Length of the first sequence (seq1): {first_key_length}")
print(f"Length of the tenth sequence (seq10): {tenth_key_length}")
print(f"Length of the last sequence ({last_key_name}): {last_key_length}")



#%%
import numpy as np 

""" Remove consecutive duplicated in the sequences"""

import numpy as np

def remove_consecutive_duplicates(s):
    # Function to remove consecutive duplicates in a string
    result = [s[0]]  # Start with the first character
    for char in s[1:]:  # Iterate over the string, starting from the second character
        if char != result[-1]:  # If the current character is not the same as the last in result
            result.append(char)  # Append the current character to the result
    return ''.join(result)

# Load the rearranged sequences
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

# Dictionary to store the shortened sequences
shortened_sequences = {}

# Iterate over each sequence, preprocess it, and store in the dictionary
for key in npz_file.files:
    seq = npz_file[key]
    if isinstance(seq, np.ndarray) and seq.shape == ():  # Convert scalar arrays to Python scalars
        seq = seq.item()
    # Apply the preprocessing function to remove consecutive duplicates
    shortened_seq = remove_consecutive_duplicates(seq)
    # Store the shortened sequence in the dictionary
    shortened_sequences[key] = shortened_seq

# Save the shortened sequences in a new npz file
output_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences.npz"
np.savez_compressed(output_path, **shortened_sequences)


#%%

""" Remove the spacer"""
import numpy as np

def remove_consecutive_duplicates_and_s(s):
    result = []  # Start with an empty list
    previous_char = None  # Keep track of the previous character
    
    for char in s:
        # Check if the current character is different from the previous character and is not 'S'
        if char != previous_char and char != 'S':
            result.append(char)  # Append the current character to the result
            previous_char = char  # Update the previous character
            
    return ''.join(result)

# Load the rearranged sequences
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

# Dictionary to store the processed sequences
shortened_sequences_without_S = {}

# Iterate over each sequence, preprocess it, and store in the dictionary
for key in npz_file.files:
    seq = npz_file[key]
    if isinstance(seq, np.ndarray) and seq.shape == ():  # Convert scalar arrays to Python scalars
        seq = seq.item()
    
    # Apply the updated preprocessing function
    shortened_seq_without_S = remove_consecutive_duplicates_and_s(seq)
    
    # Store the processed sequence in the dictionary
    shortened_sequences_without_S[key] = shortened_seq_without_S

# Save the processed sequences in a new npz file
output_path_without_S = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"
np.savez_compressed(output_path_without_S, **shortened_sequences_without_S)


#%%
import numpy as np
""" 6. Compare original vs basecalled sequence """
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']
print(original_seq)

path_rearranged = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\decoded_sequences.npz"
npz_file = np.load(path_rearranged)
# Get the keys of the npz file
keys = npz_file.files

# Print the keys
print(keys)


#%%

import numpy as np

# Load the original sequence
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['signal_train']

# Load the rearranged sequences
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

# Iterate over each sequence and print
for key in npz_file.files:
    seq = npz_file[key]
    print("Original Seq:", original_seq)
    print(f"{key}:", seq)

#%%
import numpy as np

# Load the original sequence data
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']

# Convert the numerical original_seq to a string using the provided mapping
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}  # Assuming this mapping is correct
original_seq_str = ''.join(index_to_base[base] for base in original_seq)

# Load the processed sequences
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
npz_file = np.load(path_rearranged)

# Extract seq1 and print its length
seq1 = npz_file['seq1']
if isinstance(seq1, np.ndarray) and seq1.shape == ():  # Convert scalar arrays to Python scalars
    seq1 = seq1.item()
print("Length of seq1:", len(seq1))

# Print the shape of the original sequence
print("Shape of the original sequence:", original_seq.shape)

#%%
import numpy as np

# Load the original sequence data
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
original_seq = train_data_file['rand_seq']

# Convert the numerical original_seq to a string using the provided mapping, excluding 'S'
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
original_seq_str = ''.join(index_to_base.get(base, '') for base in original_seq)  # Use .get() to handle missing keys

# Load the processed sequences that have 'S' removed and no consecutive duplicates
# Assuming you saved this new .npz file correctly after processing
path_rearranged = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\shortened_sequences_without_S.npz"
npz_file = np.load(path_rearranged)

# Extract seq1 and print its length
seq1 = npz_file['seq1']
if isinstance(seq1, np.ndarray) and seq1.shape == ():  # Convert scalar arrays to Python scalars
    seq1 = seq1.item()
print("Length of seq1 (after processing):", len(seq1))

# Print the length of the original sequence string (after removing 'S' and converting to string)
print("Length of the original sequence string:", len(original_seq_str))


#%%







""" 6. Compare original vs basecalled sequence using levenstein """
import numpy as np
from Levenshtein import distance as levenshtein_distance

# Function to decode rand_seq
def decode_rand_seq(rand_seq, index_to_base):
    return ''.join(index_to_base[base] for base in rand_seq)

# Function to compare the original sequence with each of the 32 decoded sequences
def compare_sequences(original_sequence, decoded_sequences_file):
    comparison_results = {}
    for key in decoded_sequences_file.files:
        decoded_sequence = decoded_sequences_file[key]
        if isinstance(decoded_sequence, np.ndarray) and decoded_sequence.shape == ():
            decoded_sequence = decoded_sequence.item()  # Convert to Python scalar if needed
        
        # Calculate the distance and store it
        distance = levenshtein_distance(original_sequence, decoded_sequence)
        comparison_results[key] = distance
    return comparison_results

# Load the decoded sequences
decoded_sequences_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\decoded_sequences.npz"
decoded_sequences_file = np.load(decoded_sequences_path)

# Load rand_seq and decode it
train_data_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\train_data_0.npz"
train_data_file = np.load(train_data_path)
rand_seq = train_data_file['rand_seq']
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}
original_sequence = decode_rand_seq(rand_seq, index_to_base)

# Compare the original sequence with the decoded sequences
comparison_results = compare_sequences(original_sequence, decoded_sequences_file)

# Print comparison results
for seq_key, distance in comparison_results.items():
    print(f"Edit distance for {seq_key}: {distance}")





#%%%

#older

import numpy as np
from Levenshtein import distance as levenshtein_distance

# Load the decoded sequences
decoded_sequences_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\recombined_long_read_basecalling.npz"
decoded_sequences_file = np.load(decoded_sequences_path)

# Load the original sequence (rand_seq) from train_data_0.npz
train_data_path = "C:\\Users\\manue\\MASTER_PROJECT_RNA_seq_data\\Optimize_ML_simulated_RNA_sequencing_data-main\\Optimize_ML_simulated_RNA_sequencing_data-main\\LONG_READ_training_data\\train_data_0.npz"
train_data_file = np.load(train_data_path)
rand_seq = train_data_file['signal_train']

# Convert rand_seq from numerical values to string using the same mapping
index_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'S'}
original_sequence = ''.join(index_to_base[base] for base in rand_seq)
print(original_sequence)
# Function to compare each decoded sequence with the original
def compare_sequences(decoded_sequences_file, original_sequence):
    results = {}
    for key in decoded_sequences_file.files:
        decoded_sequence = decoded_sequences_file[key]
        if isinstance(decoded_sequence, np.ndarray) and decoded_sequence.shape == ():
            decoded_sequence = decoded_sequence.item()  # Convert to Python scalar if needed
        # Calculate Levenshtein distance (edit distance) between decoded and original sequence
        distance = levenshtein_distance(decoded_sequence, original_sequence)
        results[key] = distance
    return results

# Compare and get the results
comparison_results = compare_sequences(decoded_sequences_file, original_sequence)

# Print comparison results
for seq_key, distance in comparison_results.items():
    print(f"Edit distance for {seq_key}: {distance}")






# %%
