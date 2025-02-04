import numpy as np
import vtk
import matplotlib.pyplot as plt
import os

# Specify your file path
file_directory = "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/air_wet/"

# Function to extract data from VTP file
def extract_data_from_vtp(file_path):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()

        polydata = reader.GetOutput()

        # Extract PC (assuming it's stored as point data)
        pc_array = np.array(polydata.GetPointData().GetArray("pc"))

        # Extract S_liq (assuming it's stored as point data)
        s_liq_array = np.array(polydata.GetPointData().GetArray("S_liq"))

        return pc_array, s_liq_array
    except Exception as e:
        print(f"Failed to extract data from {file_path}: {e}")

# Function to count the number of VTP files in the specified directory
def count_vtp_files(directory):
    vtp_files = [file for file in os.listdir(directory) if file.endswith(".vtp")]
    return len(vtp_files)

# Determine the number of VTP files
num_vtp_files = count_vtp_files(file_directory)

# Generate file paths dynamically
file_paths = [os.path.join(file_directory, f"2p_air_wet-{i:05d}.vtp") for i in range(num_vtp_files)]

# Specify your file name
file_name = os.path.basename(file_paths[0])[:-10]

# Extract data from each VTP file
data_arrays = [extract_data_from_vtp(file_path) for file_path in file_paths]

# Filter out None values (failed extractions)
data_arrays = [data for data in data_arrays if data is not None]

# Print the number of VTP files detected
print(f"Number of VTP files detected: {num_vtp_files}")

# Define markers and colors for each file
markers = ['o', 's', '^', 'D', 'v', '*', 'X', 'P', 'h', '+']  # Add more markers if needed
colors = plt.cm.tab10(np.linspace(0, 1, len(data_arrays)))

# Pore indices to plot
pore_indices = [4, 13, 22, 31]

# Create subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot PC against S_liq for each VTP file in subplots
for idx, pore_index in enumerate(pore_indices):
    ax = axs[idx // 2, idx % 2]  # Get subplot axis
    ax.set_xlabel('Local S$_w$ [-]', fontsize=12)
    ax.set_ylabel('Local P$_c$ [Pa]', fontsize=12)
    ax.set_title(f'Pore {pore_index}', fontsize=14)  # Use basename to get file name
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
    ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)

    for i, (file_path, (pc_array, s_liq_array)) in enumerate(zip(file_paths, data_arrays)):
        file_label = os.path.basename(file_path)[-7:-4]  # Extract the file name without extension
        marker = markers[i % len(markers)]  # Cycle through markers
        color = colors[i]  # Get color from the predefined list
        
        ax.plot(s_liq_array[pore_index], pc_array[pore_index], color=color, marker=marker, markersize=6, label=file_label)
    
    # Add legend to the current subplot
    #ax.legend(fontsize=8, loc='upper right', ncol=4)

# Add overall title to the plot
fig.suptitle(f'Comparison of Local Pore Data of ({file_name}) Across 4 Layers and {num_vtp_files-1} Timesteps', fontsize=16)

# Adjust layout
plt.tight_layout()

# Save subplots as .png file with file names
subplots_filename = os.path.join(file_directory, f"subplots_of_{file_name}.png")
plt.savefig(subplots_filename, bbox_inches='tight')
plt.close()
'''
# Save individual plots for each pore index with file names
for pore_index in pore_indices:
    # Plot PC against S_liq for each VTP file
    plt.figure(figsize=(12, 8))
    
    plt.xlabel('Local S$_w$ [-]', fontsize=14)
    plt.ylabel('Local P$_c$ [Pa]', fontsize=14)
    plt.title(f'Local P$_c$ vs. Local S$_w$ @ pore body-{pore_index} of ({file_name}) across {num_vtp_files-1} time_steps', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    for i, (file_path, (pc_array, s_liq_array)) in enumerate(zip(file_paths, data_arrays)):
        file_label = os.path.basename(file_path)[-7:-4]  # Extract the file name without extension
        marker = markers[i % len(markers)]  # Cycle through markers
        color = colors[i]  # Get color from the predefined list
        plt.plot(s_liq_array[pore_index], pc_array[pore_index], color=color, marker=marker, markersize=8, label=file_label)
    
    #plt.legend(fontsize=10, loc='upper right', ncol=4)  # Adjust ncol to control the number of columns in the legend
    
    # Save plot as .png file with file names
    plot_filename = os.path.join(file_directory, f"plot_pore_{pore_index}_of_{file_name}.png")
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory
'''