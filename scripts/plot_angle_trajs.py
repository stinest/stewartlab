
#   This script contains a function that
#   plots a line representing the
#   trajectory of a nanostar's bond angle.

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use('science')        # converts plots to LaTeX style
plt.rc("axes", labelsize=14, labelweight="bold", labelpad=10)
plt.rc("font", size=12, weight="light")   

# desired valency + condition + angle to plot
exp_condition = '0.15M'
exp_valency = '4'
angle = '12'

desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")         # extracts corresponding master file
master_dir = os.path.join(desktop_path, "ucla", "stewart_lab", "nanostar_angles")
master_file = os.path.join(master_dir, f"{exp_valency}arm_{exp_condition}_angles.csv")
master_df = pd.read_csv(master_file)
angle_data = master_df[angle]

# sim time unit = 3.06*1e-12 s
# sim time = steps*dt = 1e7*0.0005 = 5000
# sim time in units = 5000*3.06*1e-12 s = 15.3 ns
total_time_ns = 15.3
time_ns = np.linspace(0, total_time_ns, len(angle_data))

angle_columns = [col for col in master_df.columns if col != "Experiment"]       # extracts column names
colors = {'3': 'steelblue', '4': 'seagreen', '5': 'darkorange'}
color = colors.get(exp_valency)         # colors per structure
superscripts = {'3': 't', '4': 'f', '5': 'p'}
superscript = superscripts.get(exp_valency, '')     # superscripts per structure

x_label = "Simulation Time (ns)"
y_label = fr"$\theta_{{{angle}}}^{{\:{superscript}}}$ (°)"
average_angle = angle_data.mean()

plt.figure(figsize=(10, 3))
plt.plot(time_ns, angle_data, linestyle='-', linewidth=1.5, color=color, zorder=3, label=f"{exp_valency}-NS")
plt.xlabel(x_label, fontsize=18, fontweight='bold')
plt.ylabel(y_label, fontsize=18, fontweight='bold')
plt.axhline(y=average_angle, color='gray', linestyle='--', linewidth=1, label=f"{average_angle:.2f}°")
plt.xlim(0, time_ns[-1])
plt.yticks(np.arange(0, 181, 30))
plt.xticks(np.arange(0, 16, 3))
plt.gca().set_autoscale_on(False)
plt.legend(fontsize=14)

desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')       # saves plot image to desktop
filename = os.path.join(desktop_path, f'trajs_{exp_valency}.pdf')
plt.savefig(filename, bbox_inches='tight', dpi=200)