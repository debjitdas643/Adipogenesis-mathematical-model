# Distribution of cells and Fitting to gaussian distribution to find overlap region

import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture


# Get list of all files matching "output_*.out"
files = sorted(glob.glob("pulse_*.out"))

# Initialize lists to store data
cell_ids = []
cell_fates = []

# Process each file
for file in files:
    cell_id = int(file.split('_')[-1].split('.')[0])  # Extract cell number
    data = np.loadtxt(file, skiprows=1)  # Read data, skipping header

    last_expression = data[-1, 3]  # Get last PPARG expression value

    # Determine cell fate
    cell_fate = last_expression

    # Store in lists
    cell_ids.append(cell_id)
    cell_fates.append(cell_fate)

# Convert lists to NumPy arrays
cell_ids = np.array(cell_ids)
cell_fates = np.array(cell_fates)

# Save as a file
np.savetxt("actual_event.out", np.column_stack((cell_ids, cell_fates)), fmt="%d %.2f", header="Cell_ID Cell_Fate", comments="")


# Print arrays
# print("Cell IDs:", cell_ids)
# print("Cell Fates:", cell_fates)
plt.hist(cell_fates, bins=50, alpha=0.3, color='yellow', edgecolor='black')
plt.xlabel("PPARG Expression",fontsize=18)
plt.ylabel("Frequency",fontsize=18)
plt.title("Distribution of PPARG Expression",fontsize=18)
plt.show()


# Load your dataset
data = cell_fates

# Fit Gaussian Mixture Model (GMM) with 2 components
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(data.reshape(-1, 1))



# Extract mean, std, and weights
mu1, mu2 = gmm.means_.flatten()
std1, std2 = np.sqrt(gmm.covariances_).flatten()
weight1, weight2 = gmm.weights_

# Define separate x-ranges for each Gaussian
x1 = np.linspace(mu1 - 3*std1, mu1 + 3*std1, 500)
x2 = np.linspace(mu2 - 3*std2, mu2 + 3*std2, 500)

# Compute properly scaled Gaussian curves
gaussian1 = weight1 * norm.pdf(x1, mu1, std1)
gaussian2 = weight2 * norm.pdf(x2, mu2, std2)

# Find overlapping region
common_x = np.linspace(max(mu1 - 3*std1, mu2 - 3*std2), min(mu1 + 3*std1, mu2 + 3*std2), 500)
common_gaussian1 = weight1 * norm.pdf(common_x, mu1, std1)
common_gaussian2 = weight2 * norm.pdf(common_x, mu2, std2)
overlap = np.minimum(common_gaussian1, common_gaussian2)

# Find the first and last x values where overlap is significant
overlap_indices = np.where(overlap > 0.001 * max(overlap))[0]  # Threshold to ignore tiny values
first_overlap_x = common_x[overlap_indices[0]]
last_overlap_x = common_x[overlap_indices[-1]]

# Plot histogram
plt.hist(data, bins=50, density=True, alpha=0.3, color='y', edgecolor='black')

# Plot Gaussians within their relevant ranges
plt.plot(x1, gaussian1, 'r', linewidth=3, label="Gaussian 1")
plt.plot(x2, gaussian2, 'g', linewidth=3, label="Gaussian 2")

# Fill the overlapping region
plt.fill_between(common_x, overlap, color='purple', alpha=0.8, label="Overlap")

# Labels and legend
plt.title("Fitting Two Gaussians to PPARG Expression Data",fontsize=18)
plt.xlabel("PPARGa Expression",fontsize=18)
plt.ylabel("Density",fontsize=18)
plt.legend()
plt.show()

# Print overlap range
print(f"Overlap starts at x = {first_overlap_x:.1f}")
print(f"Overlap ends at x = {last_overlap_x:.1f}")

# Save data
filename = "Overlap_data.out"
np.savetxt(filename,np.column_stack((first_overlap_x,last_overlap_x)), fmt="%.2f", header="Lower_boundary Upper_boundary", comments="")
