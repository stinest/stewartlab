import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# setting files for analysis
directory = './../sims/SU24/'
top_name = directory + '1M_35C/3_5J-GUAC.top'
dat_name = directory + '1M_35C/trajectory_sim.dat'

# reading in topology+data files into data frames
df_top = pd.read_csv(top_name, delimiter=' ', names=range(4), header=0)
df_dat = pd.read_csv(dat_name, delimiter=' ', header=None, names=range(3), usecols=[0,1,2])

# VAR+FUNC DEFINITIONS
l_se = 5            # of nucleotides in sticky ends, counting unpaired base(s)
l_core = 5          # of unpaired bases at the core
N_arm = 3           # of arms in nanostar
exp = 'T = 35$^\circ$C'

N = len(df_top[0])          # of nucleotides in nanostar
l_strand = len(df_top[0][df_top[0] == 1])           # strand length
l_arm = int((l_strand - (l_se + l_core))/2)         # arm length
N_strand = int(N/l_strand)          # of strands in nanostar

N_conf = (df_dat[0].values == 't').sum()            # of configurations in the given data file
box_size = df_dat.iloc[1,2]         # box size

# function that finds core/strand indices
def find_indices(df):
    core_indice = l_arm - 1
    strand_indice = l_strand - l_se - 1
    selected_core_indices = [core_indice]
    selected_strand_indices = [strand_indice]

    for i in range(2, N_arm+1):
        for index, row in df.iterrows():
            if row[0] == i:          # enters loop at the beginning of a strand
                core_row_candidate = index + core_indice          # creates core/strand index per strand
                strand_row_candidate = index + strand_indice
                if core_row_candidate in df.index:
                    core_row = df.index.get_loc(core_row_candidate)          # gets actual indexed location
                    selected_core_indices.append(core_row)          # adds to list
                if strand_row_candidate in df.index:
                    strand_row = df.index.get_loc(strand_row_candidate)
                    selected_strand_indices.append(strand_row)
                break
        
    return selected_core_indices, selected_strand_indices
    
# function that finds center of mass by averages of core nucleotides
def calculate_COM(df, core_indices):
    # selects rows, extracts columns with coords 
    core_nucleotides = df.iloc[core_indices, -3:]
    
    # calculates average for COM coordinates
    avg_x = core_nucleotides.iloc[:, 0].mean()
    avg_y = core_nucleotides.iloc[:, 1].mean()
    avg_z = core_nucleotides.iloc[:, 2].mean()
    
    # creates/returns a data frame with calculated COM coords
    COM_coords = pd.DataFrame({'com_X': [avg_x], 'com_Y': [avg_y], 'com_Z': [avg_z]})
    return (COM_coords)

# function that calculates the bond angle between 2 strand indices
def calculate_angle(df, st1_index, st2_index):
    
    # selects row of vector nucleotide, extracts columns with coords
    st1_coords = df.iloc[strand_indices[st1_index], -3:]
    st2_coords = df.iloc[strand_indices[st2_index], -3:]
    
    # calculates COM
    COM_coords = calculate_COM(df, core_indices)
    
    # calculates vectors
    st1_vector = st1_coords - COM_coords.values[0]
    st2_vector = st2_coords - COM_coords.values[0]
    
    # converts to arrays for numpy use
    st1_array = np.array(st1_vector)
    st2_array = np.array(st2_vector)
    
    # calculates bond angle
    dot_product = np.dot(st1_array, st2_array)
    st1_magnitude = np.linalg.norm(st1_array)
    st2_magnitude = np.linalg.norm(st2_array)
    var = (dot_product / (st1_magnitude * st2_magnitude))
    bond_angle = np.degrees(np.arccos(var))
    
    return (bond_angle)

# function that calculates arm length between two indices
def calculate_armLength(df, strand_index):
    std_coords = df.iloc[strand_indices[strand_index], -3:]
    COM_coords = calculate_COM(df, core_indices)
    std_vector = std_coords - COM_coords.values[0]
    std_array = np.array(std_vector)
    return np.linalg.norm(std_array)
    
# function that makes histograms for each arm for frequency of bond angles
def make_histogram():
    graph_counter = 1      # initializes theta count for plots
    
    for i in range(N_arm-1):
        for j in range(i+1, N_arm):
            if N_arm == 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_histogram(df_calculated_angles['angles'], fr'$\theta_{graph_counter}$', fr'$P(\theta_{graph_counter})$', graph_counter)
                graph_counter += 1
            elif N_arm == 4 and i+j != 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_histogram(df_calculated_angles['angles'], fr'$\theta_{graph_counter}$', fr'$P(\theta_{graph_counter})$', graph_counter)
                graph_counter += 1   
   
# function that makes line histogram for frequency of bond angles, meant to be comparative with other plots             
def make_plot(experiment):
    calculated_angles = [calculate_angle(df, 0, 1) for df in df_conf]
    df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
    
    avg, sd = norm.fit(df_calculated_angles['angles'])
    plt.hist(df_calculated_angles['angles'], bins=5, density=True, alpha=0.6, color='white')
    
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, avg, sd)
    
    plt.plot(x, p, 'm.-', label=experiment, linewidth=1) # b+ g* m.
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$P(\theta)$')
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True, bottom=True, left=True)
    plt.minorticks_on()
    plt.legend(frameon=False)
    plt.xlim(0, 200)
    plt.ylim(0, 0.017)
    plt.show()
    
def make_histogram_armLength():
    graph_counter = 1
    
    for i in range(N_arm):
        if N_arm == 3:
            calculated_lengths = [calculate_armLength(df, i) for df in df_conf]
            df_calculated_lengths = pd.DataFrame({'lengths': calculated_lengths})
            plot_armLength(df_calculated_lengths['lengths'], f'Arm {graph_counter} Length (SU)', 'Frequency')
            plt.show()
            graph_counter += 1
        elif N_arm == 4:
            calculated_lengths = [calculate_armLength(df, i) for df in df_conf]
            df_calculated_lengths = pd.DataFrame({'lengths': calculated_lengths})
            plot_armLength(df_calculated_lengths['lengths'], f'Arm {graph_counter} Length (SU)', 'Frequency')
            plt.show()
            graph_counter += 1
    
# function that plots histograms for frequency of bond angles
def plot_histogram(data, xaxis, yaxis, counter, bins=6):
    plt.figure(figsize=(6, 5))  # Set figure size to create a square plot
    plt.hist(data, bins=bins, color='olive', edgecolor='darkgreen', alpha=0.7, align='left')
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    
    avg = np.mean(data)
    sd = np.std(data)
    plt.text(0.04, 0.96, f'$μ = {avg:.2f}$\n$σ = {sd:.2f}$', transform=plt.gca().transAxes, fontsize=15, color='black', ha='left', va='top')
    
    desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
    #filename = os.path.join(desktop_path, f'{title}_{counter}.png')
    filename = os.path.join(desktop_path, f'i_{counter}.png')
    plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    #plt.show()
    
#function that plots histograms for average arm length
def plot_armLength(data, xaxis, yaxis, bins=5):
    plt.figure(figsize=(5, 5))
    plt.hist(data, bins=bins, color='olive', edgecolor='darkgreen', alpha=0.7, align='left')
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    
    avg = np.mean(data)
    sd = np.std(data)
    plt.text(0.04, 0.96, f'$μ = {avg:.2f}$\n$σ = {sd:.2f}$', transform=plt.gca().transAxes, fontsize=15, color='black', ha='left', va='top')
    plt.show()

    
## CODE IMPLEMENTATION


# partitioning the data of each configuration
df_conf = [[] for i in range(N_conf)]
for i in range(N_conf):
    df_conf[i] = df_dat.iloc[i*N + (i+1)*3:(i+1)*N + (i+1)*3, 0:3]
    df_conf[i].reset_index(inplace=True)
    # convert dataframe to numerical values
    df_conf[i][0] = df_conf[i][0].astype(float)
    df_conf[i][1] = df_conf[i][1].astype(float)

core_indices, strand_indices = find_indices(df_top)      # finds indices

make_histogram()