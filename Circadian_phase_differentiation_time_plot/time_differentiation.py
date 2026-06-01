import numpy as np
import glob
import matplotlib.pyplot as plt

# Get list of all files matching "output_*.out"
files = sorted(glob.glob("pulse_*.out"))


threshold = int(input("Enter the threshold value obtained from ROC analysis: "))

differentiation_time = []
cell_ids = []

for file in files:
    data = np.loadtxt(file)
    time = data[:, 0]        # First column: Time
    time = time/(24*60) # time in day scale
    expression = data[:, 3]  # Forth column: PPARG expression
    cell_id = int(file.split('_')[-1].split('.')[0])  # File number as label
    crossing_event = []
    for i in range(len(expression)):
        if i != 0:
            if (threshold-expression[i])*(threshold-expression[i-1]) < 0:
                crossing_event.append(time[i])
    if len(crossing_event) != 0 and expression[-1] > (threshold+40):
        if crossing_event[-1] <= 3.0:
            differentiation_time.append(crossing_event[-1])
            cell_ids.append(cell_id)

print(len(differentiation_time))

# Save as a file
np.savetxt("Time_differentiation.out", np.column_stack((cell_ids, differentiation_time)), fmt=["%d","%.4f"])
