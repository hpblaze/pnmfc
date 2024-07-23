import numpy as np
import vtk
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET

# Specify your file paths
file_directory = "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/p1e-3ov_t3e-5/"
pvd_file_path = os.path.join(file_directory, "2p_magic_reg_het_ps1e-3ov_t3e-5.pvd")

# Function to extract timesteps from PVD file
def extract_timesteps_from_pvd(pvd_file_path):
    try:
        tree = ET.parse(pvd_file_path)
        root = tree.getroot()
        timesteps = [float(dataset.attrib["timestep"]) for dataset in root.iter("DataSet")]
        return timesteps
    except Exception as e:
        print(f"Failed to extract timesteps from {pvd_file_path}: {e}")
        return None

# Function to extract pc data from VTP file
def extract_pc_from_vtp(file_path, pore_indices):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        polydata = reader.GetOutput()
        pc_array = np.array(polydata.GetPointData().GetArray("pc"))
        return min(pc_array[pore_indices[0]], pc_array[pore_indices[1]])
    except Exception as e:
        print(f"Failed to extract pc data from {file_path}: {e}")
        return None

# Function to extract transmissibilityN data from VTP file
def extract_transmissibilityN_from_vtp(file_path, throat_index):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        polydata = reader.GetOutput()
        transmissibilityN_array = np.array(polydata.GetCellData().GetArray("transmissibilityN"))
        return transmissibilityN_array[throat_index]
    except Exception as e:
        print(f"Failed to extract transmissibilityN data from {file_path}: {e}")

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
    file_paths = [os.path.join(file_directory, f"2p_magic_reg_het_ps1e-3ov_t3e-5-{i:05d}.vtp") for i in range(num_vtp_files)]

    # Print the number of VTP files detected
    print(f"Number of VTP files detected: {num_vtp_files}")

    # Pore indices combinations to consider
    pore_index_combinations = [
        [4, 5], [13, 14], [22, 23], [31, 32], [4, 13], [13, 22], [22, 31]
    ]

    # Throat indices to plot
    throat_indices = [4, 21, 38, 55, 12, 29, 46]

    # Create subplots
    fig, axs = plt.subplots(2, 4, figsize=(18, 10))

    # Define markers and colors based on the number of timesteps
    markers = ['o', 's', '^', 'D', 'v', '*', 'X', 'P', 'h', '+']
    num_markers = len(markers)
    colors = plt.cm.tab10(np.linspace(0, 1, num_vtp_files))

    # Plot transmissibilityN against pc data for each pore index combination
    for i, (pore1, pore2) in enumerate(pore_index_combinations):
        # Calculate pc,ij values for each timestep
        pc_ij_array = [extract_pc_from_vtp(file_path, [pore1, pore2]) for file_path in file_paths]

        # Get corresponding transmissibilityN data for the throat index
        transmissibilityN_array = [extract_transmissibilityN_from_vtp(file_path, throat_indices[i]) for file_path in file_paths]

        # Plot pc,ij array against transmissibilityN data with markers based on timestep
        ax = axs[i // 4, i % 4]
        for j, (pc, transmissibility) in enumerate(zip(pc_ij_array, transmissibilityN_array)):
            ax.plot(pc, transmissibility, marker=markers[j % num_markers], color=colors[j], linestyle='', markersize=6, label=f'Timestep {j}')

        ax.set_xlabel('Local $p^{im}_{c,ij}$ [-]', fontsize=12)
        ax.set_ylabel('Local $g^{im}_{n,ij}$ [m$^3$/(Pa.s)]', fontsize=12)
        ax.set_title(f'Pore: [{pore1}, {pore2}] for throat: [{throat_indices[i]}]', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
        ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
        #ax.legend()

        # Add a red dotted line at x = 1450
        ax.axvline(x=2416.67, color='r', linestyle='--')

        # Set logarithmic scale for x-axis
        #ax.set_xscale('log')
        ax.set_xlim(0, None)

    # Adjust layout
    plt.tight_layout()

    # Save subplots as .png file with file names
    subplots_filename = os.path.join(file_directory, "gN_vs_pc_ij_subplot.png")
    plt.savefig(subplots_filename, bbox_inches='tight')

    # Show plot
    plt.show()
else:
    print("Failed to extract timesteps from PVD file. Please check the file format.")