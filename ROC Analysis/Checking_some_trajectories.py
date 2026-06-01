# Checking some trajectories

import numpy as np
import glob
import matplotlib.pyplot as plt

# Get list of all files matching "output_*.out"
files = sorted(glob.glob("pulse_*.out"))[310:320]  # First 10 files

# Plot expression vs. time for first 10 files
plt.figure(figsize=(8, 6))

for file in files:
    data = np.loadtxt(file, skiprows=1)  # Skip header row
    time = data[:, 0]        # First column: Time
    time = time/(24*60)	     # Converting to days
    expression = data[:, 3]  # Second column: PPARG expression
    plt.plot(time, expression, alpha=0.75, label=file.split('_')[-1].split('.')[0])  # File number as label

plt.xlabel("Time (Day)",fontsize=18)
plt.ylabel("PPARGa Expression",fontsize=18)
plt.title("PPARGa Expression Dynamics",fontsize=18)
plt.axis([0,4,0,2000])
plt.yticks(np.arange(0, 2201, 100))
plt.legend(title="Cell ID", loc="upper left", bbox_to_anchor=(1, 1))
plt.show()
