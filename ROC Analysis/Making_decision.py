# The distribution is density plot of last expression of each cell
# Therefore, we will get two distributions (one for preadipocyte and another for adipocyte) and there can be overlap region
# If the last expreesion is higher than upper boundary of overlap region then it belongs to positive distribution
# If the last expreesion is lower than lower boundary of overlap region then it belongs to negative distribution
# Last expression of those cells which lie in this overlap region we will check their dynamics (plot trajectory)
# If there is chance of increase in expression then it will be grouped with positive distribution otherwise negative distribution
# Decision for the cells in overlap region will be made with intuition

import numpy as np
import matplotlib.pyplot as plt
import glob

with open("Overlap_data.out","r") as f:
    next(f) # Skip header
    for line in f:
        first_overlap_x, last_overlap_x = float(line.split()[0]),float(line.split()[1])
# Get sorted list of files
files = sorted(glob.glob("pulse_*.out"))

# Step 1: Read actual fates
cell_ids = []  # List to store cell IDs
last_exp = []  # List to store last expression

with open("actual_event.out", "r") as f:
    next(f)  # Skip header
    for line in f:
        cell_id, last = int(line.split()[0]), float(line.split()[1])
        cell_ids.append(cell_id)
        last_exp.append(last)

positive_events = []
negative_events = []
for i in range(len(last_exp)):
    if last_exp[i] > last_overlap_x:
        positive_events.append(last_exp[i])
    elif last_exp[i] < first_overlap_x:
        negative_events.append(last_exp[i])
    else:
        # print(last_exp[i])
        data_to_read = files[i] #correcting indices
        data_file = np.loadtxt(data_to_read, skiprows=1)  # Skip header row
        time = data_file[:, 0]        # First column: Time
        time = time/(24*60)	      # Converting to days
        expression = data_file[:, 3]  # Second column: PPARG expression
        print(expression[-1])
        plt.plot(time, expression, label=f'Cell Id: {i+1}')  # File number as label
        plt.xlabel("Time (day)",fontsize=18)
        plt.ylabel("PPARGa Expression",fontsize=18)
        plt.axis([0,4,0,2200])
        plt.legend(loc='upper left')
        plt.show(block=False)  # Non-blocking show
        plt.pause(0.1)  # Pause to allow the plot to render

        decision = input("Is this a positive event? (y/n): ")
        plt.close()  # Close the plot after decision

        if decision.lower() == 'y':
            positive_events.append(last_exp[i])
        elif decision.lower() == 'n':
            negative_events.append(last_exp[i])
fig, ax = plt.subplots()

# Histograms
ax.hist(positive_events, bins=30, alpha=0.6, color='#D08FFF',
        edgecolor=None, label='Positive Events')
ax.hist(negative_events, bins=20, alpha=0.4, color='#00AC00',
        edgecolor=None, label='Negative Events')

# Labels
ax.set_xlabel(r'PPARG$_{\mathrm{a}}$ Expression (a.u.)',
              fontsize=18, fontname='Times New Roman', fontweight='bold')
ax.set_ylabel("Frequency",
              fontsize=18, fontname='Times New Roman', fontweight='bold')
for tick in ax.get_xticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')
    tick.set_fontsize(14)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')
    tick.set_fontsize(14)
plt.show()

# Save data
filename = "Positive_distribution.out"
np.savetxt(filename,positive_events, fmt="%.2f")

filename = "Negative_distribution.out"
np.savetxt(filename,negative_events, fmt="%.2f")
