#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

directory_path = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data'

#npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]

def plot_generated_long_read(signal_train, map_onehot, rand_seq_numeric, file_name):
    # Convert numeric sequence to string using a predefined dictionary
    base_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    rand_seq = ''.join([base_dict[num] for num in rand_seq_numeric])
    
    fig, axs = plt.subplots(3, 1, figsize=(14, 20), constrained_layout=True)

    # Plot Signal Train
    axs[0].plot(signal_train)
    axs[0].set_title('Signal Train')

    # PContour plot for map_onehot
    if map_onehot.any():
        axs[1].contourf(map_onehot.T, cmap='viridis', levels=np.linspace(0, 1, num=50))
        axs[1].set_title('Contour plot of One-hot Encoded Map')
        axs[1].set_xlabel('Time poinst')
        axs[1].set_ylabel('Nucleotide position')

    # Heatmap using Seaborn
    sns.heatmap(map_onehot.T, cmap="YlGnBu", cbar_kws={'label': 'Feature Activation'}, ax=axs[2])
    axs[2].set_title('Map Onehot Features (Bases) over time points')
    axs[2].set_xlabel('Time Points')
    axs[2].set_ylabel('Features/Bases and Spacer')
    # Adjusting the y-ticks to show all feature labels
    axs[2].set_yticks(range(map_onehot.shape[1]))
    axs[2].set_yticklabels([f'Feature {i}' for i in range(map_onehot.shape[1])]) #change to letter, spacer later

    plt.suptitle(file_name)
    plt.show()
    # sequence as a string
    print(f"Random Sequence for {file_name}: {rand_seq}")
    print(f"Signal Train for {file_name}:\n", signal_train)
    print(f"Map Onehot for {file_name}:\n", map_onehot)

# Loop through each file and plot (currently one 1 npz file)
npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]
for npz_file in npz_files:
    file_path = os.path.join(directory_path, npz_file)
    data = np.load(file_path)
    signal_train = data['signal_train'][0] 
    map_onehot = data['map_onehot'][0]  
    rand_seq_numeric = data['rand_seq']  

plot_generated_long_read(signal_train, map_onehot, rand_seq_numeric, npz_file)


# %%
