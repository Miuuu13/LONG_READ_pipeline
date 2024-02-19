#%%
import numpy as np
import matplotlib.pyplot as plt
import os

# directory containing my .npz files
directory_path = r'C:\Users\manue\MASTER_PROJECT_RNA_seq_data\Optimize_ML_simulated_RNA_sequencing_data-main\Optimize_ML_simulated_RNA_sequencing_data-main\LONG_READ_training_data\Simulated_LONG_READ_data_FEB24\Simulated_data_n2_w800_b32'

# List all .npz files in this directory
npz_files = [f for f in os.listdir(directory_path) if f.endswith('.npz')]

""" Plotting function """
def plot_data(signal_train, map_onehot, file_name):
    fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axs[0].plot(signal_train, label='Signal Train')
    axs[0].set_title('Signal Train')
    axs[0].legend()

    for i in range(map_onehot.shape[1]):
        axs[1].plot(map_onehot[:, i], label=f'Channel {i}')
    axs[1].set_title('One-hot Encoded Map')
    axs[1].legend()

    plt.suptitle(file_name)
    plt.xlabel('Time Points')
    plt.show()

# Loop through each file and plot
for npz_file in npz_files:
    file_path = os.path.join(directory_path, npz_file)
    data = np.load(file_path)
    signal_train = data['signal_train'][0]  # Assuming you want to plot the first sequence in each file
    map_onehot = data['map_onehot'][0]  # Adjust indices as needed

    plot_data(signal_train, map_onehot, npz_file)
