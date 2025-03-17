
#   This script contains the foundational analysis to generate the necessary 
#   data frames from a given nanostar simulation and defines functions to 
#   calculate bond angles for use in other scripts.

import pandas as pd
import numpy as np
import os       # for master file of all sims

# user-defined variables + files
salt = 1
temp = 37
N_sims = 3          # of repeated simulations
l_se = 7            # of nucleotides in sticky ends, counting unpaired base(s)
l_core = 2          # of unpaired bases at the core
directory = './../sims/DNA_3m/GCTAGC/2bp/'
experiment = f"{salt}M_{temp}C"        # sim parameters/folder name
top_name = f"{directory}{experiment}/3m1_6NT1_2bp.top"

# reading in topology + trajectory files into data frames
df_top = pd.read_csv(top_name, delimiter=' ', names=range(4), header=0) 
df_dat_list = []        # reads in trajectory files for all repeated sims
for i in range(1, N_sims+1):
    dat_names = f"{directory}{experiment}/trajectory_sim_{i}.dat"
    df_dats = pd.read_csv(dat_names, delimiter=' ', header=None, names=range(3), usecols=[0, 1, 2])
    df_dat_list.append(df_dats)
filtered_dfs = []       # averages numeric rows (uses energy from first sim)
for df in df_dat_list:
    numeric_rows = ~df[0].astype(str).str.match(r'^(t|b|E|=)$')         # filters numeric rows
    numeric_df = df[numeric_rows].astype(float)         # converts to float for averaging
    filtered_dfs.append(numeric_df)
concat_df = pd.concat(filtered_dfs, axis=0)         # concatenates filtered numeric rows 
averaged_numeric_df = concat_df.groupby(concat_df.index).mean()         # averages
non_numeric_rows = ~df_dat_list[0][0].astype(str).str.match(r'^(t|b|E|=)$')
non_numeric_df = df_dat_list[0][~non_numeric_rows]      # preserves non-numeric rows
df_dat = pd.concat([non_numeric_df, averaged_numeric_df]).sort_index()      # combines preserved + averaged rows

# calculated parameters
N = len(df_top[0])          # of nucleotides in nanostar
l_strand = len(df_top[0][df_top[0] == 1])           # strand length
l_arm = int((l_strand - (l_se + l_core))/2)         # arm length
N_arm = int(N/l_strand)          # of arms in nanostar
N_conf = (df_dat[0].values == 't').sum()            # of configurations in the given data file
label = fr"S = {salt}M Na$^+$, T = {temp}°C"

# partitioning the data of each configuration
df_conf = [[] for i in range(N_conf)]
for i in range(N_conf):
    df_conf[i] = df_dat.iloc[i*N + (i+1)*3:(i+1)*N + (i+1)*3, 0:3]
    df_conf[i].reset_index(inplace=True)
    # convert dataframe to numerical values
    df_conf[i][0] = df_conf[i][0].astype(float)
    df_conf[i][1] = df_conf[i][1].astype(float)         # returns df_conf with averaged df_dats values
    
        
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

core_indices, strand_indices = find_indices(df_top)      # finds indices
df_calculated_angles = pd.DataFrame()       # initializes an empty df to store angles
allowed_pairs_5arm = [(0, 1), (0, 2), (1, 4), (2, 3), (3, 4)]
custom_labels_4arm = ["14", "12", "34", "23"]       # predefined labels for 4-armed according to .top
custom_labels_5arm = ["15", "12", "45", "23", "34"]     # predefined labels for 5-armed according to .top
label_index = 0     # tracks custom labels

for i in range(N_arm - 1):
    for j in range(i + 1, N_arm):
        if N_arm == 3:
            column_label = f"{i + 1}{j + 1}"        # since labels are accurate to strand labels for 3-armed
            calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
            df_calculated_angles[column_label] = calculated_angles
        elif N_arm == 4 and i+j != 3:
            column_label = custom_labels_4arm[label_index]
            calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
            df_calculated_angles[column_label] = calculated_angles
            label_index += 1
        elif N_arm == 5:
            if (i, j) in allowed_pairs_5arm:            # filters to right angles
                column_label = custom_labels_5arm[label_index]
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles[column_label] = calculated_angles
                label_index += 1

# call to save angle data in a master file
def save_to_file():
    # adding experiment's bond angles to a master file
    desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")         
    master_dir = os.path.join(desktop_path, "ucla", "stewart_lab", "nanostar_angles")        # path to master files
    os.makedirs(master_dir, exist_ok=True)      # confirms directory exists

    df_calculated_angles["Experiment"] = experiment         # adds new column for new experiment
    salt_file = os.path.join(master_dir, f"DNA_{N_arm}arm_{salt}M_angles.csv")         # new master files for new parameters
    temp_file = os.path.join(master_dir, f"DNA_{N_arm}arm_{temp}C_angles.csv")
    
    append_to_master(salt_file, df_calculated_angles)
    append_to_master(temp_file, df_calculated_angles)
# appends to existing or creates new master file per experiment
def append_to_master(path, data):
    if os.path.exists(path):
        master_df = pd.read_csv(path)       # loads existing master file
        if data["Experiment"].iloc[0] in master_df["Experiment"].values:        # makes sure same runs are not being duplicated
            print(f"experiment '{data['Experiment'].iloc[0]}' already exists in {path}. skipping append.")
            return
        master_df = pd.concat([master_df, data], ignore_index=True)     # appends new experiment to file
    else:
        master_df = data        # creates new master file if none exist
    master_df.to_csv(path, index=False)         # saves file
    print(f"saved to: {path}")      # confirms save

# accessors
def get_experiment():
    return experiment
def get_label():
    return label
def get_angles():
    return df_calculated_angles

save_to_file()