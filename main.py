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
file_path = r"C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_basecalled\train_data_0_fragment_1.npz"
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
