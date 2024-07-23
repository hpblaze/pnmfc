import numpy as np
import vtk
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET

def extract_timesteps_from_pvd(pvd_file_path):
    """Extract timesteps from a PVD file."""
    try:
        tree = ET.parse(pvd_file_path)
        root = tree.getroot()
        timesteps = [float(dataset.attrib["timestep"]) for dataset in root.iter("DataSet")]
        return timesteps
    except Exception as e:
        print(f"Failed to extract timesteps from {pvd_file_path}: {e}")
        return None

def extract_transmissibilityW_from_vtp(file_path):
    """Extract transmissibilityW data from a VTP file."""
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        polydata = reader.GetOutput()
        transmissibilityW_array = np.array(polydata.GetCellData().GetArray("transmissibilityW"))
        return transmissibilityW_array
    except Exception as e:
        print(f"Failed to extract transmissibilityW data from {file_path}: {e}")
        return None

def extract_transmissibilityN_from_vtp(file_path):
    """Extract transmissibilityN data from a VTP file."""
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        polydata = reader.GetOutput()
        transmissibilityN_array = np.array(polydata.GetCellData().GetArray("transmissibilityN"))
        return transmissibilityN_array
    except Exception as e:
        print(f"Failed to extract transmissibilityN data from {file_path}: {e}")
        return None

def count_vtp_files(directory, base_filename):
    """Count the number of VTP files in the specified directory."""
    return len([file for file in os.listdir(directory) if file.startswith(base_filename) and file.endswith(".vtp")])

# Specify your file paths
file_directory = "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_sink/time_1sec"
base_filename = "2p_im_magic_time_1sec"
pvd_file_path = os.path.join(file_directory, f"{base_filename}.pvd")

# Extract timesteps from PVD file
timesteps = extract_timesteps_from_pvd(pvd_file_path)

if timesteps:
    # Determine the number of VTP files
    num_vtp_files = count_vtp_files(file_directory, base_filename)

    # Generate file paths dynamically
    file_paths = [os.path.join(file_directory, f"{base_filename}-{i:05d}.vtp") for i in range(num_vtp_files)]

    # Extract transmissibilityW and transmissibilityN data from each VTP file
    transmissibilityW_data_arrays = [extract_transmissibilityW_from_vtp(file_path) for file_path in file_paths]
    transmissibilityN_data_arrays = [extract_transmissibilityN_from_vtp(file_path) for file_path in file_paths]

    # Filter out None values (failed extractions)
    transmissibilityW_data_arrays = [data for data in transmissibilityW_data_arrays if data is not None]
    transmissibilityN_data_arrays = [data for data in transmissibilityN_data_arrays if data is not None]

    # Print the number of VTP files detected
    print(f"Number of VTP files detected: {len(transmissibilityW_data_arrays)}")

    # Define colors, markers, and marker size
    color_w = 'blue'
    color_n = 'red'
    marker_w = 'o'
    marker_n = 's'
    marker_size = 5  # Adjust this value to set the marker size

    # Throat indices to plot
    throat_indices = [12, 29, 46]

    # Calculate the number of rows and columns dynamically
    num_plots = len(throat_indices)
    num_cols = 3
    num_rows = 1  # num_plots // num_cols + (1 if num_plots % num_cols != 0 else 0)

    # Create subplots for all plots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 5))

    # Ensure axs is always a 2D array
    if num_rows == 1 and num_cols == 1:
        axs = np.array([axs])
    elif num_rows == 1 or num_cols == 1:
        axs = np.ravel(axs)

    # Plot transmissibilityW and transmissibilityN against physical time for each throat index
    for idx, throat_index in enumerate(throat_indices):
        ax = axs[idx]  # Get subplot axis
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Local $g^{im}_{w,ij}$ [m$^3$/(Pa.s)]', fontsize=12, color=color_w)
        ax.set_title(f'Throat {throat_index}', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.tick_params(axis='x', which='major', direction='out', length=6, width=1, labelsize=10)
        ax.tick_params(axis='y', which='major', direction='out', length=6, width=1, labelsize=10, labelcolor=color_w)
        ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
        
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)

        # Create a twin y-axis to plot transmissibilityN
        ax2 = ax.twinx()
        ax2.set_ylabel('Local $g^{im}_{n,ij}$ [m$^3$/(Pa.s)]', fontsize=12, color=color_n)
        ax2.tick_params(axis='y', labelcolor=color_n)

        w_values = []
        n_values = []
        times = []

        for i, (transmissibilityW_array, transmissibilityN_array) in enumerate(zip(transmissibilityW_data_arrays, transmissibilityN_data_arrays)):
            time_step = timesteps[i]  # Physical time corresponding to each timestep
            times.append(time_step)
            w_values.append(transmissibilityW_array[throat_index])
            n_values.append(transmissibilityN_array[throat_index])

        ax.plot(times, w_values, color=color_w, marker=marker_w, markersize=marker_size, linestyle='None', label=f'Throat {throat_index} (W)')
        ax2.plot(times, n_values, color=color_n, marker=marker_n, markersize=marker_size, linestyle='None', label=f'Throat {throat_index} (N)')

    # Adjust layout
    plt.tight_layout()

    # Save subplots as .png file with file names
    subplots_filename = os.path.join(file_directory, "merged_gW_and_gN_vs_realtime.png")
    plt.savefig(subplots_filename, bbox_inches='tight')
    plt.show()
else:
    print("Failed to extract timesteps from PVD file. Please check the file format.")
