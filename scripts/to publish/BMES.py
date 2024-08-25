import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

# color scheme: firebrick, royalblue, yellowgreen
import scienceplots
plt.style.use('science')

# setting files for analysis
directory = './../sims/S24/'
top_name = directory + '0.5M_35C/3arm_4SE.top'
dat_name = directory + '0.5M_35C/trajectory_sim.dat'
title = 'nanostar = 3arm_4SE(-GUAC)'
experiment = 'T = 35\u00B0C'
experiment_number = 3

# reading in topology+data files into data frames
df_top = pd.read_csv(top_name, delimiter=' ', names=range(4), header=0)
df_dat = pd.read_csv(dat_name, delimiter=' ', header=None, names=range(3), usecols=[0,1,2])

# VAR+FUNC DEFINITIONS
l_se = 5            # of nucleotides in sticky ends, counting unpaired base(s)
l_core = 2          # of unpaired bases at the core
N_arm = 3           # of arms in nanostar

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
    
    
    
## FUNCTIONS
    
    
    
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
def plot_histogram(data, xaxis, yaxis, counter, bins=5):
    plt.figure(figsize=(5, 5))
    plt.hist(data, bins=bins, color='royalblue', edgecolor='black', alpha=0.8, align='left', zorder=3)
    plt.xlabel(xaxis, fontsize=12, fontweight='bold')
    plt.ylabel(yaxis, fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.grid(True, linestyle='--', alpha=0.3, zorder=1)
    
    avg = np.mean(data)
    sd = np.std(data)
    plt.text(0.04, 0.96, f'$\\mu = {avg:.2f}$\n$\\sigma = {sd:.2f}$', transform=plt.gca().transAxes, fontsize=15, color='black', ha='left', va='top')
    
    desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
    filename = os.path.join(desktop_path, f'histogram_{counter}.pdf')
    plt.savefig(filename, bbox_inches='tight', dpi=200)
  
# function that makes plots for bond angle across time
def make_noise():
    graph_counter = 1      # initializes theta count for plots
    
    for i in range(N_arm-1):
        for j in range(i+1, N_arm):
            if N_arm == 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_noise(df_calculated_angles['angles'], 'configurations', fr'$\theta_{graph_counter}$', graph_counter)
                graph_counter += 1
            elif N_arm == 4 and i+j != 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_noise(df_calculated_angles['angles'], 'configurations', fr'$\theta_{graph_counter}$', graph_counter)
                graph_counter += 1  
def plot_noise(data, xaxis, yaxis, counter):
    plt.figure(figsize=(8, 3))
    plt.plot(data, linestyle='-', linewidth=1, color='firebrick', zorder=3)
    plt.xlim(0, 1000)
    plt.ylim(0, max(data) * 1.1)
    plt.xlabel(xaxis, fontsize=12, fontweight='bold')
    plt.ylabel(yaxis, fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.grid(True, linestyle='--', alpha=0.3, zorder=1)
        
    x = np.arange(len(data))
    coeffs = np.polyfit(x, data, 1)
    line_of_best_fit = np.polyval(coeffs, x)
    plt.plot(x, line_of_best_fit, linestyle='--', color='black', zorder=3)
    
    desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
    filename = os.path.join(desktop_path, f'noise_{counter}.pdf')
    plt.savefig(filename, bbox_inches='tight', dpi=200)

# function that makes heat maps for bond angle across all sims
def make_3d():
    df_theta1 = pd.DataFrame(data_theta1)
    df_theta2 = pd.DataFrame(data_theta2)
    df_theta3 = pd.DataFrame(data_theta3)
    
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    save_path = os.path.join(desktop_path, '3d_plots')
    os.makedirs(save_path, exist_ok=True)
    angles = [(0, 30), (90, 30), (0, 90)]
    
    #plot_3d(df_theta1, r'$\theta_1$', angles, save_path, 'theta_1')
    #plot_3d(df_theta2, r'$\theta_2$', angles, save_path, 'theta_2')
    #plot_3d(df_theta3, r'$\theta_3$', angles, save_path, 'theta_3')
    
    plot_combined_3d([df_theta1, df_theta2, df_theta3], ['firebrick', 'royalblue', 'yellowgreen'], angles, save_path, 'combined_thetas')
def plot_3d(df, zaxis, angles, save_path, prefix):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = df['x']
    y = df['y']
    z = df['z']
    X, Y = np.meshgrid(np.linspace(x.min(), x.max(), 30), np.linspace(y.min(), y.max(), 30))
    Z = griddata((x, y), z, (X, Y), method='cubic')

    # Plot the surface
    ax.plot_surface(X, Y, Z, cmap="autumn_r", edgecolor='k', lw=0.5, rstride=1, cstride=1)
    
    # Add contours
    ax.contour(X, Y, Z, 10, lw=3, cmap="autumn_r", linestyles="solid", offset=-1)
    ax.contour(X, Y, Z, 10, lw=3, colors="k", linestyles="solid")

    ax.set_xlabel('temperature (C)')
    ax.set_ylabel('salt concentration (M)')
    ax.set_zlabel(zaxis)
    ax.set_xticks([10, 25, 35])
    ax.set_yticks([0.1, 0.5, 1])
    
    desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
    filename = os.path.join(desktop_path, f'{prefix}.pdf')
    plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    for angle, elevation in angles:
        ax.view_init(elevation, angle)
        filename = os.path.join(save_path, f'{prefix}_angle_{angle}_elev_{elevation}.pdf')
        plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    plt.close(fig)
def plot_combined_3d(dfs, colors, angles, save_path, prefix):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111, projection='3d')

    for df, color in zip(dfs, colors):
        x = df['x']
        y = df['y']
        z = df['z']
        X, Y = np.meshgrid(np.linspace(x.min(), x.max(), 30), np.linspace(y.min(), y.max(), 30))
        Z = griddata((x, y), z, (X, Y), method='cubic')
        
        ax.plot_surface(X, Y, Z, color=color, alpha=0.5, edgecolor='k', lw=0.5, rstride=1, cstride=1, zorder=3)

    ax.set_xlabel('Temp (°C)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Salt Conc (M)', fontsize=12, fontweight='bold')
    ax.set_zlabel('Bond Angle (°)', fontsize=12, fontweight='bold')
    ax.set_xticks([10, 25, 35])
    ax.set_yticks([0.1, 0.5, 1])
    
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.grid(True, linestyle='--', alpha=0.3, zorder=1)
    
    filename = os.path.join(save_path, f'{prefix}.pdf')
    plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    for angle, elevation in angles:
        ax.view_init(elevation, angle)
        filename = os.path.join(save_path, f'{prefix}_angle_{angle}_elev_{elevation}.pdf')
        plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    plt.close(fig)

# function that makes line histogram for bond angle frequency, meant for overlay          
def make_hist_trend():
    graph_counter = 1      # initializes theta count for plots
    
    for i in range(N_arm-1):
        for j in range(i+1, N_arm):
            if N_arm == 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_hist_trend(df_calculated_angles['angles'], fr'$\theta_{graph_counter}$', fr'$P(\theta_{graph_counter})$', graph_counter)
                graph_counter += 1
            elif N_arm == 4 and i+j != 3:
                calculated_angles = [calculate_angle(df, i, j) for df in df_conf]
                df_calculated_angles = pd.DataFrame({'angles': calculated_angles}) 
                plot_hist_trend(df_calculated_angles['angles'], fr'$\theta_{graph_counter}$', fr'$P(\theta_{graph_counter})$', graph_counter)
                graph_counter += 1   
def plot_hist_trend(df, xaxis, yaxis, counter):
    plt.figure(figsize=(8, 5))
    hist_counts, bin_edges = np.histogram(df, bins=20)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    hist_frequencies = hist_counts / len(df)
    
    interpolator = interp1d(bin_centers, hist_frequencies, kind='cubic', fill_value="extrapolate")      # fills in more data
    x_dense = np.linspace(bin_edges[0], bin_edges[-1], 40)
    y_dense = interpolator(x_dense)
    
    if experiment_number == 1:
        plotstyle = 'bs-'
    elif experiment_number == 2:
        plotstyle = 'gD-'
    elif experiment_number == 3:
        plotstyle = 'mo-'
        
    plt.plot(x_dense, y_dense, plotstyle, markerfacecolor='none', zorder=4, label=experiment)

    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.ylim(0, 0.16)
    plt.xlim(20, 180)
    plt.legend()
    
    desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
    filename = os.path.join(desktop_path, f'histtrend_{counter}.pdf')
    plt.savefig(filename, bbox_inches='tight', dpi=300)


    
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

data_theta1 = {
    'x': [10, 10, 10, 25, 25, 25, 35, 35, 35],
    'y': [0.1, 0.5, 1, 0.1, 0.5, 1, 0.1, 0.5, 1],
    'z': [117.9, 111.04, 120.66, 112.04, 117.04, 134.43, 113.32, 118.91, 115.8] }
data_theta2 = {
    'x': [10, 10, 10, 25, 25, 25, 35, 35, 35],
    'y': [0.1, 0.5, 1, 0.1, 0.5, 1, 0.1, 0.5, 1],
    'z': [112.92, 99.61, 97.94, 113.63, 120.12, 94.87, 116.47, 106.32, 107.29] }
data_theta3 = {
    'x': [10, 10, 10, 25, 25, 25, 35, 35, 35],
    'y': [0.1, 0.5, 1, 0.1, 0.5, 1, 0.1, 0.5, 1],
    'z': [114.41, 132.05, 124.14, 111.79, 103.5, 114.7, 113.69, 115.6, 119.65] }

#make_histogram()
#make_noise()
#make_3d()
#make_hist_trend()
