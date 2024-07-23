# -*- coding: utf-8 -*-
"""Pc_sw curve plot from vtp.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1JP6iOMmhsyMekuu45fn1xAX7M3eyM1hI
"""

import numpy as np
import vtk
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET

# Specify your file paths
file_directory = "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/3300nm_throat_len/"
pvd_file_path = os.path.join(file_directory, "2p_small_len_throat_3300nm.pvd")

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

# Function to extract S_liq data from VTP file
def extract_s_liq_from_vtp(file_path):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()

        polydata = reader.GetOutput()

        # Extract S_liq (assuming it's stored as point data)
        s_liq_array = np.array(polydata.GetPointData().GetArray("S_liq"))

        return s_liq_array
    except Exception as e:
        print(f"Failed to extract S_liq data from {file_path}: {e}")
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
    file_paths = [os.path.join(file_directory, f"2p_small_len_throat_3300nm-{i:05d}.vtp") for i in range(num_vtp_files)]

    # Extract S_liq data from each VTP file
    s_liq_data_arrays = [extract_s_liq_from_vtp(file_path) for file_path in file_paths]

    # Filter out None values (failed extractions)
    s_liq_data_arrays = [s_liq_data for s_liq_data in s_liq_data_arrays if s_liq_data is not None]

    # Print the number of VTP files detected
    print(f"Number of VTP files detected: {num_vtp_files}")

    # Define markers and colors for each file
    markers = ['o', 's', '^', 'D', 'v', '*', 'X', 'P', 'h', '+']  # Add more markers if needed
    colors = plt.cm.tab10(np.linspace(0, 1, len(s_liq_data_arrays)))

    # Pore indices to plot
    pore_indices = [4, 13, 22, 31, 40, 49, 58, 67, 76]

    # Calculate the number of rows and columns dynamically
    num_plots = len(pore_indices)
    num_cols = 5
    num_rows = (num_plots + num_cols - 1) // num_cols

    # Create subplots for all plots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 4 * num_rows))

    # Plot S_liq against timestep for each VTP file in subplots
    for idx, pore_index in enumerate(pore_indices):
        ax = axs[idx // num_cols, idx % num_cols]  # Get subplot axis
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Local S$_w$ [-]', fontsize=12)
        ax.set_title(f'Pore {pore_index}', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
        ax.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)

        for i, s_liq_array in enumerate(s_liq_data_arrays):
            marker = markers[i % len(markers)]  # Cycle through markers
            color = colors[i]  # Get color from the predefined list
            ax.plot(timesteps[i], s_liq_array[pore_index], color=color, marker=marker, markersize=6, label=f'Time Step {i}')

    # Remove empty subplots
    for idx in range(num_plots, num_rows * num_cols):
        fig.delaxes(axs.flatten()[idx])

    # Adjust layout
    plt.tight_layout()

    # Save subplots as .png file with file names
    subplots_filename = os.path.join(file_directory, "subplots_S_liq_with_realtime.png")
    plt.savefig(subplots_filename, bbox_inches='tight')
    plt.show()
else:
    print("Failed to extract timesteps from PVD file. Please check the file format.")
