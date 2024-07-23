import numpy as np
import vtk
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET

# Specify your file paths
file_directory = "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_1e-04/"
pvd_file_path = os.path.join(file_directory, "2p_throat_rad_1e-04.pvd")

# Function to extract timesteps from PVD file
def extract_timesteps_from_pvd(pvd_file_path):
    try:
        tree = ET.parse(pvd_file_path)
        root = tree.getroot()
        timesteps = []

        for dataset in root.iter("DataSet"):
            timestep = float(dataset.attrib["timestep"])
            timesteps.append(timestep)

        return timesteps
    except Exception as e:
        print(f"Failed to extract timesteps from {pvd_file_path}: {e}")
        return None

# Function to extract p_gas data from VTP file
def extract_p_gas_from_vtp(file_path):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()

        polydata = reader.GetOutput()

        # Extract p_gas (assuming it's stored as point data)
        p_gas_array = np.array(polydata.GetPointData().GetArray("p_gas"))

        return p_gas_array
    except Exception as e:
        print(f"Failed to extract p_gas data from {file_path}: {e}")
        return None
        
# Function to count the number of VTP files in the specified directory
def count_vtp_files(directory):
    vtp_files = [file for file in os.listdir(directory) if file.endswith(".vtp")]
    return len(vtp_files)

# Determine the number of VTP files
num_vtp_files = count_vtp_files(file_directory)

# Extract timesteps from PVD file
timesteps = extract_timesteps_from_pvd(pvd_file_path)

if timesteps:
    # Generate file paths dynamically
    file_paths = [os.path.join(file_directory, f"2p_throat_rad_1e-04-{i:05d}.vtp") for i in range(num_vtp_files)]

    # Extract p_gas data from each VTP file
    p_gas_data_arrays = [extract_p_gas_from_vtp(file_path) for file_path in file_paths]

    # Filter out None values (failed extractions)
    p_gas_data_arrays = [p_gas_data for p_gas_data in p_gas_data_arrays if p_gas_data is not None]

    # Print the number of VTP files detected
    print(f"Number of VTP files detected: {num_vtp_files}")

    # Define markers and colors for each file
    markers = ['o', 's', '^', 'D', 'v', '*', 'X', 'P', 'h', '+']  # Add more markers if needed
    colors = plt.cm.tab10(np.linspace(0, 1, len(p_gas_data_arrays)))

    # Pore indices to plot
    pore_indices = [4, 13, 22, 31]

    # Create subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # Plot p_gas against timestep for each VTP file in subplots
    for idx, pore_index in enumerate(pore_indices):
        ax = axs[idx // 2, idx % 2]  # Get subplot axis
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Local $p_{n}$ [Pa]', fontsize=12)
        ax.set_title(f'Pore {pore_index}', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
        ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)
        #ax.ylim(99999,100001)

        for i, p_gas_array in enumerate(p_gas_data_arrays):
            marker = markers[i % len(markers)]  # Cycle through markers
            color = colors[i]  # Get color from the predefined list
            ax.plot(timesteps[i], p_gas_array[pore_index], color=color, marker=marker, markersize=6, label=f'Time Step {i}')

    # Adjust layout
    plt.tight_layout()

    # Save subplots as .png file with file names
    subplots_filename = os.path.join(file_directory, "subplots_p_gasvsRealtime.png")
    plt.savefig(subplots_filename, bbox_inches='tight')
    plt.show()
else:
    print("Failed to extract timesteps from PVD file. Please check the file format.")
'''
# Save individual plots for each pore index with file names
for pore_index in pore_indices:
    # Plot p_gas against time_step for each VTP file
    plt.figure(figsize=(12, 8))
    
    plt.xlabel('Time step', fontsize=14)
    plt.ylabel('Local P$_n$ [-]', fontsize=14)
    plt.title(f'Local P$_n$ vs. timestep @ pore body-{pore_index}', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    for i, p_gas_array in enumerate(p_gas_data_arrays):
        time_step = i + 1  # Time step starts from 0
        marker = markers[i % len(markers)]  # Cycle through markers
        color = colors[i]  # Get color from the predefined list
        plt.plot(time_step, p_gas_array[pore_index], color=color, marker=marker, markersize=8, label=f'Time Step {i}')
    
    plt.legend(fontsize=10, loc='upper left', ncol = 4)  # Adjust legend location and font size
    
    # Save plot as .png file with file names
    plot_filename = os.path.join(file_directory, f"plot_pore_{pore_index}_p_gas_with_legend.png")
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()  # Close the current figure to avoid overlapping plots
'''