#%%

# not cmpleted
import os
import numpy as np 
import matplotlib.pyplot as plt

base_path = os.getcwd()

def compare_sequences(generated_seq, basecalled_seq):
    """
    Compare two sequences and identify mismatches.

    Parameters:
    - generated_seq: The original sequence generated.
    - basecalled_seq: The sequence obtained from basecalling.

    Returns:
    - mismatches: A list of tuples, each containing the index and the mismatching characters (generated, basecalled).
    """
    def compare_sequences(generated_seq, basecalled_seq):
        mismatches = []
    # Ensure both sequences are 1D arrays; try to flatten them if they're not
        if generated_seq.ndim > 1:
            generated_seq = generated_seq.flatten()
        if basecalled_seq.ndim > 1:
            basecalled_seq = basecalled_seq.flatten()
    
    # Check if lengths match, adjust if necessary (this might need a more sophisticated approach depending on your needs)
        min_length = min(len(generated_seq), len(basecalled_seq))
        generated_seq = generated_seq[:min_length]
        basecalled_seq = basecalled_seq[:min_length]

        for i in range(min_length):
            gen_char = generated_seq[i]
            call_char = basecalled_seq[i]
            if gen_char != call_char:
                mismatches.append((i, gen_char, call_char))
        return mismatches

def plot_mismatches(generated_seq, basecalled_seq, mismatches):
    """
    Plot the original and basecalled sequences, highlighting mismatches in red.

    Parameters:
    - generated_seq: The original generated sequence.
    - basecalled_seq: The basecalled sequence.
    - mismatches: A list of mismatches as returned by compare_sequences.
    """
    plt.figure(figsize=(20, 3))
    for i, char in enumerate(generated_seq):
        color = 'red' if i in [m[0] for m in mismatches] else 'black'
        plt.text(i, 1, char, color=color, ha='center')
    
    for i, char in enumerate(basecalled_seq):
        color = 'red' if i in [m[0] for m in mismatches] else 'black'
        plt.text(i, 0, char, color=color, ha='center')

    plt.axis('off')
    plt.show()


original_seq_path = os.path.join(base_path, 'LONG_READ_training_data')
npz_files_original = [f for f in os.listdir(original_seq_path) if f.endswith('.npz')]
if npz_files_original:
    npz_file_path = os.path.join(original_seq_path, npz_files_original[0])  # Load the first file
    data = np.load(npz_file_path)
    original_seq = data['rand_seq']
else:
    print("No .npz files found in the specified basecalled directory.")

basecalled_seq_path = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
basecalled_npz_files = [f for f in os.listdir(basecalled_seq_path) if f.endswith('.npz')]
if basecalled_npz_files:
    basecalled_npz_file_path = os.path.join(basecalled_seq_path, basecalled_npz_files[0])  # Load the first file
    basecalled_data = np.load(basecalled_npz_file_path)
    basecalled_seq = basecalled_data['basecalling']
else:
    print("No .npz files found in the specified basecalled directory.")

mismatches = compare_sequences(original_seq, basecalled_seq)
plot_mismatches(original_seq, basecalled_seq, mismatches)


#%%























import os
import numpy as np
import matplotlib as plt

base_path = os.getcwd()

def compare_sequences(generated_seq, basecalled_seq):
    """
    Compare two sequences and identify mismatches.

    Parameters:
    - generated_seq: The original sequence generated.
    - basecalled_seq: The sequence obtained from basecalling.

    Returns:
    - mismatches: A list of tuples, each containing the index and the mismatching characters (generated, basecalled).
    """
    mismatches = []
    for i, (gen_char, call_char) in enumerate(zip(generated_seq, basecalled_seq)):
        if gen_char != call_char:
            mismatches.append((i, gen_char, call_char))
    return mismatches


def plot_mismatches(generated_seq, basecalled_seq, mismatches):
    """
    Plot the original and basecalled sequences, highlighting mismatches in red.

    Parameters:
    - generated_seq: The original generated sequence.
    - basecalled_seq: The basecalled sequence.
    - mismatches: A list of mismatches as returned by compare_sequences.
    """
    plt.figure(figsize=(20, 3))
    for i, char in enumerate(generated_seq):
        color = 'red' if i in [m[0] for m in mismatches] else 'black'
        plt.text(i, 1, char, color=color, ha='center')
    
    for i, char in enumerate(basecalled_seq):
        color = 'red' if i in [m[0] for m in mismatches] else 'black'
        plt.text(i, 0, char, color=color, ha='center')

    plt.axis('off')
    plt.show()


original_seq_path = os.path.join(base_path, 'LONG_READ_training_data')
npz_files_original = [f for f in os.listdir(original_seq_path) if f.endswith('.npz')]
if npz_files_original:
    npz_file_path = os.path.join(original_seq_path, npz_files_original[0])  # Load the first file
    data = np.load(npz_file_path)
    original_seq = data['rand_seq']
else:
    print("No .npz files found in the specified basecalled directory.")

basecalled_seq_path = os.path.join(base_path, 'LONG_READ_training_data_basecalled')
basecalled_npz_files = [f for f in os.listdir(basecalled_seq_path) if f.endswith('.npz')]
if basecalled_npz_files:
    basecalled_npz_file_path = os.path.join(basecalled_seq_path, basecalled_npz_files[0])  # Load the first file
    basecalled_data = np.load(basecalled_npz_file_path)
    basecalled_seq = basecalled_data['basecalling']
else:
    print("No .npz files found in the specified basecalled directory.")

mismatches = compare_sequences(original_seq, basecalled_seq)
plot_mismatches(original_seq, basecalled_seq, mismatches)






#%%
import matplotlib.pyplot as plt

def compare_sequences(original_seq, basecalled_seq):
    mismatches = []
    for i, (o, c) in enumerate(zip(original_seq, basecalled_seq)):
        if o != c:
            mismatches.append((i, o, c))
    return mismatches



def plot_sequence_comparison(original_seq, basecalled_seq, mismatches):
    # Setup plot
    plt.figure(figsize=(20, 5))
    ax = plt.gca()
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Convert sequences to strings if they are in numeric format
    base_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    original_seq_str = ''.join([base_dict[b] if b in base_dict else b for b in original_seq])
    called_seq_str = ''.join([base_dict[b] if b in base_dict else b for b in basecalled_seq])

    # Plotting
    for i, (o, c) in enumerate(zip(original_seq_str, called_seq_str)):
        color = 'red' if (i, o, c) in mismatches else 'black'
        ax.text(i, 0.7, o, ha='center', va='center', color=color)
        ax.text(i, 0.3, c, ha='center', va='center', color=color)

    plt.xlim(0, len(original_seq))
    plt.show()






#%%
def compare_sequences(original_seq, predicted_seq):
    """
    Compare the original and predicted sequences.
    
    Args:
    - original_seq: The original sequence (numpy array).
    - predicted_seq: The predicted sequence (numpy array).
    
    Returns:
    - similarity_score: A simple similarity score or metric.
    """
    # Assuming original_seq and predicted_seq are numpy arrays of the same shape
    # and contain sequence information encoded in a comparable format.
    
    # A simple comparison could be the fraction of matching positions:
    match_count = np.sum(original_seq == predicted_seq)
    total_count = original_seq.size
    
    similarity_score = match_count / total_count
    return similarity_score
