import numpy as np
import matplotlib.pyplot as plt
import os

# Ensure npz_files and directory_path are defined as before

def plot_data(signal_train, map_onehot, rand_seq_numeric, file_name):
    # Convert numeric sequence to string
    base_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    rand_seq = ''.join([base_dict[num] for num in rand_seq_numeric])
    
    fig, axs = plt.subplots(3, 1, figsize=(12, 12), constrained_layout=True)

    # Plot Signal Train
    axs[0].plot(signal_train)
    axs[0].set_title('Signal Train')

    # Plot Filled Contour for map_onehot
    if map_onehot.any():
        axs[1].contourf(map_onehot.T, cmap='viridis', levels=np.linspace(0, 1, num=50))
        axs[1].set_title('Filled Contour of One-hot Encoded Map')
        axs[1].set_xlabel('Time Point')
        axs[1].set_ylabel('Nucleotide Position')

    # Plot Heatmap for map_onehot
    axs[2].imshow(map_onehot.T, aspect='auto', cmap='hot')
    axs[2].set_title('Heatmap of One-hot Encoded Map')
    axs[2].set_xlabel('Time Point')
    axs[2].set_ylabel('Nucleotide Position')

    plt.suptitle(file_name)
    plt.show()

    # Print the sequence as a string
    print(f"Random Sequence for {file_name}: {rand_seq}")

# Loop through each file and plot, including the new changes
npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]
for npz_file in npz_files:
    file_path = os.path.join(directory_path, npz_file)
    data = np.load(file_path)
    signal_train = data['signal_train'][0]  # Assuming you want to plot the first sequence in each file
    map_onehot = data['map_onehot'][0]  # Adjust indices as needed
    rand_seq_numeric = data['rand_seq']  # Assuming rand_seq is saved as numeric values

    plot_data(signal_train, map_onehot, rand_seq_numeric, npz_file)

