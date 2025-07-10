# removes extra confs

import os   

desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
input_path = os.path.join(desktop_path, "DNA_trajectory_sim_3.dat")
output_path = os.path.join(desktop_path, "trajectory_sim_3.dat")

def filter_dna_trajectory(input_file, output_file, step_interval=100000):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        write_block = False
        for line in infile:
            if line.startswith("t = "):
                timestep = int(line.split('=')[1].strip())
                write_block = (timestep % step_interval == 0)
            if write_block:
                outfile.write(line)

filter_dna_trajectory(input_path, output_path)

