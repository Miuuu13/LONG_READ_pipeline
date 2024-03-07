#%%
import tensorflow as tf
import numpy as np
import os

def predict_and_save_basecalling(input_dir, output_dir, model_path):
    # Load the trained model
    model = tf.keras.models.load_model(model_path)
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through each .npz file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.npz'):
            file_path = os.path.join(input_dir, file_name)
            
            # Load the fragment data
            with np.load(file_path) as data:
                signal_train = data['signal_train']  # Input signal
                map_onehot = data['map_onehot']  # Mapping information
                rand_seq = data['rand_seq']  # Original sequence
                
            # Predict basecalling
            predictions = model.predict(signal_train, batch_size=32)
            
            # Save the prediction and original data into a new .npz file
            output_file_path = os.path.join(output_dir, f"basecalled_{file_name}")
            np.savez(output_file_path, 
                     signal_train=signal_train, 
                     map_onehot=map_onehot, 
                     rand_seq=rand_seq, 
                     basecalling=predictions)
            # Predictions are saved as numerical representations to be decoded later

# Example usage
base_path = os.getcwd()

input_dir = os.path.join(base_path, 'LONG_READ_training_data_splitted')
output_dir = os.path.join(base_path, 'LONG_READ_training_data_basecalled')

model_path = os.path.join(base_path, 'UNET_LSTM_n2_w1200_64f_e20_b32_22FEB24.h5')

predict_and_save_basecalling(input_dir, output_dir, model_path)



#%%
"""
import tensorflow as tf
import numpy as np
import os

def predict_and_save_basecalling(input_dir, output_dir, model_path):
    model = tf.keras.models.load_model(model_path)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through each .npz file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.npz'):
            file_path = os.path.join(input_dir, file_name)
            
            # Load the .npz file
            with np.load(file_path) as data:
                x_data = data['signal_train']  # Adjust the key if necessary
                
            # Predict basecalling
            predictions = model.predict(x_data, batch_size=32)
        
            
            # Save the prediction and original data into a new .npz file
            output_file_path = os.path.join(output_dir, f"basecalled_{file_name}")
            np.savez(output_file_path, basecalling=predictions, original_data=x_data)

# directories
base_path = os.getcwd()

input_dir = os.path.join(base_path, r'LONG_READ_training_data_splitted')
output_dir = os.path.join(base_path, r'LONG_READ_training_data_basecalled')

model_path = os.path.join(base_path, r'UNET_LSTM_n2_w1200_64f_e20_b32_22FEB24.h5')

"""
# %%
