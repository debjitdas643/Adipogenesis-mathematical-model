import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks
from scipy.stats import gaussian_kde
import matplotlib as mpl

def get_phase_at_time(query_time, t_array, phase_array):
    if query_time < t_array[0] or query_time > t_array[-1]:
        return None  # Out of bounds
    return np.interp(query_time, t_array, phase_array)

def list_appending(main_list, append_list):
    for i in append_list:
        main_list.append(i)
    return main_list

# Get list of all files matching "output_*.out"
files = sorted(glob.glob("pulse_*.out"))
#print(files[0])

time_file = np.loadtxt("Time_differentiation.out")
differentiation_time = time_file[:, 1]
cell_ids_list = time_file[:, 0]
cell_id = []
circadian_phase = []


for i in cell_ids_list:
    cell_id.append(int(i))
# the above steps are so that cell_ids can be treated as integer
#print(cell_id)

selected_files = []
#files[i] for i in cell_id
for i in cell_id:
#    print(i)
#    print(files[i-1])
    selected_files.append(files[i-1]) # by adjusting the index
    
for i in range(len(selected_files)):
    data = np.loadtxt(selected_files[i])
    time = data[:, 0]        # First column: Time
    time = time/(24*60) # time in day scale
    rev = data[:, 15] # RevErb

    smoothed = savgol_filter(rev, window_length=51, polyorder=2)
    
#    plt.plot(time,smoothed)
#    plt.show()
    # Find maxima and minima
    peak_idxs, _ = find_peaks(smoothed, prominence=300)
    trough_idxs, _ = find_peaks(-smoothed, prominence=300)

    # Now convert to actual time + value pairs
    global_peaks = [(time[j], smoothed[j]) for j in peak_idxs]
    global_troughs = [(time[j], smoothed[j]) for j in trough_idxs]

    # Extract just the time values from peaks and troughs
    peak_times = np.array([t for t, _ in global_peaks])
    trough_times = np.array([t for t, _ in global_troughs])

    # Label them as 'peak' or 'trough'
    extrema_times = np.concatenate((peak_times, trough_times))
    extrema_types = np.array(['peak'] * len(peak_times) + ['trough'] * len(trough_times))

    # Sort them by time
    sorted_indices = np.argsort(extrema_times)
    extrema_times = extrema_times[sorted_indices]
    extrema_types = extrema_types[sorted_indices]

    # Initialize the phase and time array
    all_phase = []
    all_time = []

    # Interpolate phase between each peak-trough or trough-peak pair
    for k in range(len(extrema_times) - 1):
        t_start = extrema_times[k]
        t_end = extrema_times[k + 1]
        kind_start = extrema_types[k]
        kind_end = extrema_types[k + 1]

        t_segment = [t for t in time if t_start <= t <= t_end]
#        print(t_segment)
#        print(k)
#        print(t_start,kind_start,t_end,kind_end)

        phase_segment = []
        if kind_start == 'peak' and kind_end == 'trough':
        # Falling phase: 0 to pi
            phase_at_start = 0.0
            phase_at_end = np.pi
        elif kind_start == 'trough' and kind_end == 'peak':
        # Rising phase: pi to 2pi
            phase_at_start = np.pi
            phase_at_end = 2*np.pi
        
        for t in t_segment:
            phase = phase_at_start + ((phase_at_end - phase_at_start)/(t_end - t_start))*(t - t_start) # Linear Interpolation
            phase_segment.append(phase)

        all_time = list_appending(all_time,t_segment)
        all_phase = list_appending(all_phase,phase_segment)
#    print(all_time)
#    print(len(all_time))
#    print(i)

    phase_at_diff = get_phase_at_time(differentiation_time[i], all_time, all_phase)
    circadian_phase.append(phase_at_diff)
#    print(cell_id[i], differentiation_time[i], phase_at_diff)


# Save as a file
np.savetxt("Diff_time_phase.out", np.column_stack((cell_id, differentiation_time, circadian_phase)), fmt=["%d","%.4f", "%.4f"])

# Forcing all math + regular text to Times
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.weight'] = 'bold'

# phase plot
plt.hist(circadian_phase, bins=40, color='pink', alpha=0.75, edgecolor='black')
plt.xlabel("Circadian Phase (Rev)", fontsize=18, fontweight='bold')
plt.xticks(
    ticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],  # Positions at these phase values
    labels=[r'$\mathrm{0}$', r'$\mathrm{\frac{\pi}{2}}$', 
        r'$\mathrm{\pi}$', r'$\mathrm{\frac{3\pi}{2}}$', 
        r'$\mathrm{2\pi}$'],  # Labels in LaTeX format
fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Cell Count", fontsize=18, fontweight='bold')
plt.title("Phase plot", fontsize=20, fontweight='bold')
plt.tight_layout()
plt.savefig("phase_dist.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()


# diff time plot
plt.hist(differentiation_time, bins=40, color='pink', alpha=0.75, edgecolor='black')
plt.xlabel("Differentiation time (Day)", fontsize=18, fontweight='bold')
plt.ylabel("Cell Count", fontsize=18, fontweight='bold')
plt.xlim([0,3.05])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Time plot", fontsize=20, fontweight='bold')
plt.tight_layout()
plt.savefig("time_dist.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()

kde = gaussian_kde(differentiation_time)
xi = np.linspace(0, 3, 25)   # domain similar to xlim
f = kde(xi)*752
plt.fill_between(xi, f, linewidth=2, color='pink', alpha=0.45)
plt.xlabel("Differentiation time (Day)", fontsize=18, fontweight='bold')
plt.ylabel("Cell Count", fontsize=18, fontweight='bold')
plt.xlim([0,3.05])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Time plot", fontsize=20, fontweight='bold')
plt.tight_layout()
plt.savefig("time_smooth_dist.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()



# phase time plot
differentiation_time = np.array(differentiation_time)
rev_phase = np.array(circadian_phase)

xy = np.vstack([differentiation_time, rev_phase])
z = gaussian_kde(xy)(xy)

idx = z.argsort()
differentiation_time, rev_phase, z = differentiation_time[idx], rev_phase[idx], z[idx]

plt.scatter(differentiation_time, rev_phase, c=z, s=10, cmap="turbo", edgecolor=None, alpha=0.8)
# Colorbar
cbar = plt.colorbar()
cbar.set_label("Density", fontsize=18, fontweight='bold')

# Make colorbar tick labels match
cbar.ax.tick_params(labelsize=14)


plt.xlabel("Differentiation time (Day)", fontsize=18, fontweight='bold')
plt.xlim([0,3.05])
plt.yticks(
    ticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],  # Positions at these phase values
    labels=[r'$\mathrm{0}$', r'$\mathrm{\frac{\pi}{2}}$', 
        r'$\mathrm{\pi}$', r'$\mathrm{\frac{3\pi}{2}}$', 
        r'$\mathrm{2\pi}$'],  # Labels in LaTeX format
fontsize=14)
plt.xticks(fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.ylabel("Circadian Phase (Rev)", fontsize=18, fontweight='bold')
plt.title("Phase-Time plot", fontsize=20, fontweight='bold')
plt.tight_layout()
plt.savefig("phase_time_plot.svg", format="svg", bbox_inches="tight")
plt.show()
plt.close()

