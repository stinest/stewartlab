# Stewart Lab
> This repository contains the analysis codes + files I wrote/used in generating data for the Stewart Lab!

## Python Scripts

#### `scripts/code checks/check_angle`
This script contains the code that confirms the calculation for bond angles, ensuring the average remains around 360~ degrees.
#### `scripts/code checks/check_sim`
This script contains the code that confirms simulation trajectories, ensuring that the 
movement between two nucleotides within a single configuration (dr) does not exceed that of those nucleotides between consecutive configurations (d1 and d2).
#### `scripts/W24_analysis.py`
This script contains the analysis code for:
- **Average Bond Angles:** generates histograms to analyze average angle distributions
- **Average Arm Lengths:** generates histograms to analyze average arm length distributions
#### `scripts/SU24_analysis.py`
This script contains the analysis code for:
- **Bond Angles Over Time:** generates noise plots to analyze average bond angles over time (simulation progression)
- **Simulated Bond Angles:** generates 3D heat maps to analyze average angle distributions across all simulated temperatures/salt concentrations
