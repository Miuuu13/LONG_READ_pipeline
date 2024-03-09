#%%import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# directory containing my .npz files
directory_path = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data_splitted'

# List all .npz files in this directory
npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]

# Assuming npz_files and directory_path are defined as before

def plot_data(signal_train, map_onehot, rand_seq_numeric, file_name):
    # Convert numeric sequence to string using a predefined dictionary
    base_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    rand_seq = ''.join([base_dict[num] for num in rand_seq_numeric])
    
    # Start plotting
    fig, axs = plt.subplots(3, 1, figsize=(14, 20), constrained_layout=True)

    # Plot Signal Train
    axs[0].plot(signal_train)
    axs[0].set_title('Signal Train')

    # Plot Filled Contour for map_onehot
    if map_onehot.any():
        axs[1].contourf(map_onehot.T, cmap='viridis', levels=np.linspace(0, 1, num=50))
        axs[1].set_title('Filled Contour of One-hot Encoded Map')
        axs[1].set_xlabel('Time Point')
        axs[1].set_ylabel('Nucleotide Position')

    # Plot Heatmap for map_onehot using Seaborn
    sns.heatmap(map_onehot.T, cmap="YlGnBu", cbar_kws={'label': 'Feature Activation'}, ax=axs[2])
    axs[2].set_title('Map Onehot Features Over Time')
    axs[2].set_xlabel('Time Points')
    axs[2].set_ylabel('Features')
    # Adjusting the y-ticks to show all feature labels
    axs[2].set_yticks(range(map_onehot.shape[1]))
    axs[2].set_yticklabels([f'Feature {i}' for i in range(map_onehot.shape[1])])

    plt.suptitle(file_name)
    plt.show()

    # Print the sequence as a string
    print(f"Random Sequence for {file_name}: {rand_seq}")

# Loop through each file and plot
npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]
for npz_file in npz_files:
    file_path = os.path.join(directory_path, npz_file)
    data = np.load(file_path)
    signal_train = data['signal_train'][0]  
    map_onehot = data['map_onehot'][0]  
    rand_seq_numeric = data['rand_seq']  

    plot_data(signal_train, map_onehot, rand_seq_numeric, npz_file)

# %%
