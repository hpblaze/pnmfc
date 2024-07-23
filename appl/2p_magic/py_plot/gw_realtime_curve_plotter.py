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

# Function to extract transmissibilityW data from VTP file
def extract_transmissibilityW_from_vtp(file_path):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()

        polydata = reader.GetOutput()

        # Extract transmissibilityW (assuming it's stored as cell data)
        transmissibilityW_array = np.array(polydata.GetCellData().GetArray("transmissibilityW"))

        return transmissibilityW_array
    except Exception as e:
        print(f"Failed to extract transmissibilityW data from {file_path}: {e}")

# Function to count the number of VTP files in the specified directory
def count_vtp_files(directory):
    vtp_files = [file for file in os.listdir(directory) if file.endswith(".vtp")]
    return len(vtp_files)

# Extract timesteps from PVD file
timesteps = extract_timesteps_from_pvd(pvd_file_path)

if timesteps:
    # Determine the number of VTP files
    num_vtp_files = count_vtp_files(file_directory)

    # Generate file paths dynamically
    file_paths = [os.path.join(file_directory, f"2p_throat_rad_1e-04-{i:05d}.vtp") for i in range(num_vtp_files)]

    # Extract transmissibilityW data from each VTP file
    transmissibilityW_data_arrays = [extract_transmissibilityW_from_vtp(file_path) for file_path in file_paths]

    # Filter out None values (failed extractions)
    transmissibilityW_data_arrays = [data for data in transmissibilityW_data_arrays if data is not None]

    # Print the number of VTP files detected
    print(f"Number of VTP files detected: {len(transmissibilityW_data_arrays)}")

    # Define markers and colors for each file
    markers = ['o', 's', '^', 'D', 'v', '*', 'X', 'P', 'h', '+']  # Add more markers if needed
    colors = plt.cm.tab10(np.linspace(0, 1, len(transmissibilityW_data_arrays)))

    # Throat indices to plot
    throat_indices = [4, 21, 38, 55, 12, 29, 46]

    # Calculate the number of rows and columns dynamically
    num_plots = len(throat_indices)
    num_cols = 4
    num_rows = num_plots // num_cols + (1 if num_plots % num_cols != 0 else 0)

    # Create subplots for all plots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 4*num_rows))

    # Plot transmissibilityW against physical time for each throat index
    for idx, throat_index in enumerate(throat_indices):
        ax_row = idx // num_cols
        ax_col = idx % num_cols
        ax = axs[ax_row, ax_col]  # Get subplot axis
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Local $g^{im}_{w,ij}$ [m$^3$/(Pa.s)]', fontsize=12)
        ax.set_title(f'Throat {throat_index}', fontsize=14)  # Use basename to get file name
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
        ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)

        for i, transmissibilityW_array in enumerate(transmissibilityW_data_arrays):
            time_step = timesteps[i]  # Physical time corresponding to each timestep
            color = colors[i % len(colors)]  # Get color from the predefined list

            ax.plot(time_step, transmissibilityW_array[throat_index], color=color, marker=markers[i % len(markers)], markersize=6, label=f'Throat {throat_index}')

    # Adjust layout
    plt.tight_layout()

    # Save subplots as .png file with file names
    subplots_filename = os.path.join(file_directory, "transmissibilityW_vs.realtime.png")
    plt.savefig(subplots_filename, bbox_inches='tight')
    plt.show()
else:
    print("Failed to extract timesteps from PVD file. Please check the file format.")
