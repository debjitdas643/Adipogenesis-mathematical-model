# ROC Plot and Threshold determination

import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

# Forcing all math + regular text to Times
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.weight'] = 'bold'

positive_events = []
negative_events = []

with open("Positive_distribution.out","r") as f:
    for line in f:
        expression = float(line.split()[0])
        positive_events.append(expression)
        
with open("Negative_distribution.out","r") as f:
    for line in f:
        expression = float(line.split()[0])
        negative_events.append(expression)
        
thresholds = list(range(0,2401,1))

TP_list, TN_list, FP_list, FN_list = [], [], [], []

for threshold in thresholds:
    TP = TN = FP = FN = 0
    for i in positive_events:
        if i > threshold:
            TP += 1
        else:
            FN += 1
    for i in negative_events:
        if i < threshold:
            TN += 1
        else:
            FP += 1

    TP_list.append(TP)
    TN_list.append(TN)
    FP_list.append(FP)
    FN_list.append(FN)

TPR = []
FPR = []
for i in range(len(TP_list)):
    TPR.append(TP_list[i] / (TP_list[i] + FN_list[i]))
    FPR.append(FP_list[i] / (FP_list[i] + TN_list[i]))



roc_auc = auc(FPR, TPR)
print(f"AUC: {roc_auc:.4f}")


distance_list = []

threshold_ideal = [0.0,1.0]

for i in range(len(thresholds)):
    distance = np.sqrt((FPR[i]-threshold_ideal[0])**2+(TPR[i]-threshold_ideal[1])**2)
    distance_list.append(distance)
threshold_index = np.argmin(distance_list)

optimum_threshold = thresholds[threshold_index]
print(optimum_threshold)


# Create scatter plot
fig, ax = plt.subplots()

# ROC scatter points
ax.scatter(FPR, TPR, c=thresholds, cmap='turbo', edgecolors='none', s=12)

# --- Add ONE threshold point ---
x0 = FPR[threshold_index]     #FPR value at threshold
y0 = TPR[threshold_index]     #TPR value at threshold
ax.scatter(x0, y0, color='magenta', alpha=0.75, s=80, edgecolors='black', zorder=5, label='Threshold Point')
# ---------------------------------------------------------

# Colorbar
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('Threshold level in PPARG', fontsize=18, fontweight='bold', fontname='Times New Roman')

cbar.ax.tick_params(labelsize=14)   # size
for tick in cbar.ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

# Labels
ax.set_xlabel('FPR', fontsize=18, fontweight='bold', fontname='Times New Roman')
ax.set_ylabel('TPR', fontsize=18, fontweight='bold', fontname='Times New Roman')
#ax.set_title("ROC Plot", fontsize=18, fontweight='bold', fontname='Times New Roman')

# Tick formatting
for tick in ax.get_xticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')
    tick.set_fontname('Times New Roman')

for tick in ax.get_yticklabels():
    tick.set_fontsize(14)
    tick.set_fontweight('bold')
    tick.set_fontname('Times New Roman')

# Optional legend
ax.legend(fontsize=14, frameon=False)

plt.savefig("roc_2d.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()


# Create a figure
fig = plt.figure()

# Add 3D axes to the figure
ax = fig.add_subplot(111, projection='3d')

# Create a 3D scatter plot
ax.scatter(FPR, TPR, thresholds, c=thresholds, cmap='turbo', edgecolors='none', s=8, marker='o') 

ax.scatter(x0, y0, optimum_threshold, color='magenta', alpha=0.85, s=60, edgecolors='black', depthshade=False, label='Threshold Point')

# Draw a 2D overlay point (always on top)
ax.plot([x0], [y0], [optimum_threshold],
        marker='o',
        markersize=8, alpha=0.8,
        markerfacecolor='magenta',
        markeredgecolor='black',
        linestyle='None',
        zorder=10000)   # force to top

# Set labels for axes
ax.set_xlabel('FPR', fontsize=14, fontweight='bold', fontname='Times New Roman')
ax.set_ylabel('TPR', fontsize=14, fontweight='bold', fontname='Times New Roman')
ax.set_zlabel('Threshold', fontsize=14, fontweight='bold', fontname='Times New Roman')

# Tick formatting
for tick in ax.get_xticklabels():
    tick.set_fontsize(12)
    tick.set_fontweight('bold')
    tick.set_fontname('Times New Roman')

for tick in ax.get_yticklabels():
    tick.set_fontsize(12)
    tick.set_fontweight('bold')
    tick.set_fontname('Times New Roman')

for tick in ax.get_zticklabels():
    tick.set_fontsize(12)
    tick.set_fontweight('bold')
    tick.set_fontname('Times New Roman')

ax.legend(fontsize=13, loc='upper right', frameon=False)

# Show the plot
plt.savefig("roc_3d.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()

