#%%
import numpy as np
from scipy import signal
import os
import pandas as pd

""" Script to generate LONG READS (6000 time points), noise level 2 %, batch_size = 32 """

base_path = os.getcwd() # working directory # in main
path_save_long_read = os.path.join(base_path, r'LONG_READ_training_data')
# Check if the save path exists, if not, create it
if not os.path.exists(path_save_long_read ):
    os.makedirs(path_save_long_read )

# Acess kmer table to get values for the 5-mers
kmer_info = pd.read_csv('template_median68pA_200mv.txt', delim_whitespace=True)
kmer_info = kmer_info.values

def Generate_sim_signal(N_batch, batch_size , path, kmer_info, seq_len , time_points):
    """ Function to generate the LONG READ """
    
    seq_len = seq_len + 5
    labels = 5
    kmer_length = 5

    for j in range (N_batch):
        print(j)
        Sim_raw_signal_tot = np.zeros([batch_size,time_points])
        Map_raw_signal_tot = np.zeros([batch_size,time_points,labels])

        for k in range(batch_size):
            # generate random time step to read k-mer
            Rand_seq = np.random.randint(0,4, seq_len)
            rand_seq_base =np.empty(seq_len, dtype=str)
            Sim_raw_signal = np.zeros([time_points])
            Map_raw_signal_S = np.zeros([time_points,labels])

            # dictionary version one
            base_dict = { 0:"A", 1:"C", 2:"G", 3:"T"}
            base_dict_2 = { "A":0, "C":1, "G":2, "T":3}

            for i in range(len(Rand_seq)):
                rand_seq_base[i] = base_dict[Rand_seq[i]]
            probe_raw_data = 0
            mu, sigma = 3.686, 0.4254

            for i in range(seq_len - kmer_length + 1):
                #obtain k-mer Raw-signal value
                if i > 0:
                    #find the k-mer value in the table and associate the specific signal value
                    Single_kmer_vector = rand_seq_base[i: i + kmer_length]
                    Single_kmer = ''.join(rand_seq_base[i: i + kmer_length])
                    probe_kmer = np.where(kmer_info[:,0] == Single_kmer)[0]

                    # k-mer signal depends on a gaussian distrbution center on a specific value
                    # generate here the k-mer distribution and create the k-mer signal from it
                    kmer_center = kmer_info[probe_kmer[0],1]
                    kmer_std = (1/4)*kmer_info[probe_kmer[0],2]
                    Signal_kmer = np.random.normal(kmer_center, kmer_std)

                    # from the distribution, obtain the number of time points that the k-mer is present in the signal
                    N_step = int(np.random.lognormal(mu, sigma))

                    while N_step < 20:
                        N_step = int(np.random.lognormal(mu, sigma))

                    Sim_raw_signal[probe_raw_data : probe_raw_data + N_step] = Signal_kmer

                    Spc = 3
                    if i == 1:

                        if np.random.randint(2) == 0:
                            Map_raw_signal_S[probe_raw_data : probe_raw_data + N_step - Spc, base_dict_2[Single_kmer_vector[-1]]] = 1
                            Map_raw_signal_S[probe_raw_data + N_step - Spc : probe_raw_data + N_step, -1] = 1

                        else:
                            Map_raw_signal_S[probe_raw_data + Spc : probe_raw_data + N_step -Spc, base_dict_2[Single_kmer_vector[-1]]] = 1
                            #this is for adding blank signal
                            Map_raw_signal_S[0 : Spc, -1] = 1
                            Map_raw_signal_S[probe_raw_data + N_step -Spc : probe_raw_data + N_step, -1] = 1

                    else:

                        Map_raw_signal_S[probe_raw_data + Spc : probe_raw_data + N_step -Spc, base_dict_2[Single_kmer_vector[-1]]] = 1

                        #this is for adding blank signal
                        Map_raw_signal_S[probe_raw_data : probe_raw_data + Spc, -1] = 1
                        Map_raw_signal_S[probe_raw_data + N_step -Spc : probe_raw_data + N_step, -1] = 1

                    probe_raw_data += N_step

            # filter frequency
            filtered_signal = signal.savgol_filter(Sim_raw_signal, 11, 1) #was 10, 1 must be odd

            # add noise, 2 % for LONG READ
            noise_level = 0.02 
            Sim_signal_final = filtered_signal + np.max(filtered_signal)*np.random.normal(0,noise_level, len(filtered_signal))

            # code for the median. for dorado they use median
            median = np.median(Sim_signal_final)
            std_med = np.median(np.abs(Sim_signal_final-median))*1.4826 + np.finfo(np.float32).eps
            Sim_signal_final = (Sim_signal_final - median)/std_med

            Sim_raw_signal_tot[k] = Sim_signal_final #filtered_signal
            Map_raw_signal_tot[k] = Map_raw_signal_S

        file_name = "train_data_{}.npz".format(j);

        np.savez_compressed(os.path.join(path,file_name), signal_train=Sim_raw_signal_tot, map_onehot=Map_raw_signal_tot, rand_seq = Rand_seq)



# %%
