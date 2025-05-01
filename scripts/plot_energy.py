
#   This script contains a function that
#   plots a line representing the
#   energy trajectory of a simulation.

#   Written for studying multi-NS systems.

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scienceplots

plt.style.use('science')        # converts plots to LaTeX style
plt.rc("axes", labelsize=14, labelweight="bold", labelpad=10)
plt.rc("font", size=12, weight="light")   

valency = '4'
connections = '-1'
sim = 'mc'

directory = '~/Documents/oxRNA/sims'
experiment = f"multi-RNA_{valency}m"
structure = f"{valency}m_6NT_2bp"
energy_file = f"{directory}/{experiment}/{structure}/{valency}m{connections}/energy-{sim}.dat"

colors = {'mc': 'yellowgreen', 'md': 'lightblue', 'sim': 'orange'}
color = colors.get(sim)
axes = {'mc': 'Steps', 'md': 'Time (SU)', 'sim': 'Time (SU)'}
axis = axes.get(sim)

energy_df = pd.read_csv(energy_file, delim_whitespace=True, header=None, usecols=[0, 1])
energy_df.columns = ['Time', 'PE']
avg_pe = energy_df['PE'].mean()

plt.figure(figsize=(6, 4))
plt.plot(energy_df['Time'], energy_df['PE'], linewidth=1.5, color=color)
plt.xlabel(axis)
plt.ylabel('Energy (SU)')
plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

plt.axhline(avg_pe, color='crimson', linestyle='--', linewidth=0.8, label=fr'$\mu = {avg_pe:.2f}$')
plt.legend(
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0.,
    frameon=False
)

desktop_path = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
filename = os.path.join(desktop_path, f'{valency}m{connections}_energy-{sim}.pdf')
plt.savefig(filename, bbox_inches='tight', dpi=200)
