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

def count_vtp_files(directory, base_filename):
    """Count the number of VTP files in the specified directory."""
    return len([file for file in os.listdir(directory) if file.startswith(base_filename) and file.endswith(".vtp")])

def extract_data_from_dataset(directory, pvd_file, base_filename):
    """Extract data from a given dataset."""
    pvd_file_path = os.path.join(directory, pvd_file)
    timesteps = extract_timesteps_from_pvd(pvd_file_path)
    if not timesteps:
        return None, None
    
    num_vtp_files = count_vtp_files(directory, base_filename)
    file_paths = [os.path.join(directory, f"{base_filename}-{i:05d}.vtp") for i in range(num_vtp_files)]
    
    transmissibilityW_data_arrays = [extract_transmissibilityW_from_vtp(file_path) for file_path in file_paths]
    transmissibilityW_data_arrays = [data for data in transmissibilityW_data_arrays if data is not None]
    
    return timesteps, transmissibilityW_data_arrays

# Directories and file patterns for each dataset
datasets = [
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_1e-04/", "pvd": "2p_throat_rad_1e-04.pvd", "base_filename": "2p_throat_rad_1e-04"},
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_5e-04/", "pvd": "2p_throat_rad_5e-04.pvd", "base_filename": "2p_throat_rad_5e-04"},
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_8e-05/", "pvd": "2p_throat_rad_8e-05.pvd", "base_filename": "2p_throat_rad_8e-05"},
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_8e-05/", "pvd": "2p_throat_rad_8e-05.pvd", "base_filename": "2p_throat_rad_8e-05"}
]

# Initialize lists to hold the merged data
merged_timesteps = []
merged_transmissibilityW_data = []
labels = []

# Loop through each dataset and extract data
for dataset in datasets:
    timesteps, transmissibilityW_data_arrays = extract_data_from_dataset(
        dataset["directory"], dataset["pvd"], dataset["base_filename"])
    
    if timesteps and transmissibilityW_data_arrays:
        merged_timesteps.append(timesteps)
        merged_transmissibilityW_data.append(transmissibilityW_data_arrays)
        labels.append(dataset["base_filename"])
    else:
        print(f"Failed to process dataset: {dataset['pvd']}")

# Plotting the merged data for throat 46
throat_index = 46
plt.figure(figsize=(10, 6))
plt.xlabel('Time [s]', fontsize=12)
plt.ylabel('Local $g^{im}_{w,ij}$ [m$^3$/(Pa.s)]', fontsize=12)
plt.title(f'Throat {throat_index} Comparison', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
plt.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
plt.gca().xaxis.set_tick_params(width=1.5)
plt.gca().yaxis.set_tick_params(width=1.5)


#{"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/2p_throat_rad_8e-05/", "pvd": "2p_throat_rad_8e-05.pvd", "base_filename": "2p_throat_rad_8e-05"},
# Use different colors and markers for each dataset , 'black' , '*'
colors = ['blue', 'green', 'red','black']
markers = ['o', 's', '^', '*']

for i, (timesteps, transmissibilityW_data_arrays) in enumerate(zip(merged_timesteps, merged_transmissibilityW_data)):
    for j, transmissibilityW_array in enumerate(transmissibilityW_data_arrays):
        time_step = timesteps[j]
        plt.plot(time_step, transmissibilityW_array[throat_index], color=colors[i], marker=markers[i], markersize=3, label=f'{labels[i]}' if j == 0 else "")

plt.legend()
plt.tight_layout()

# Save the plot as a .png file
plot_filename = os.path.join(datasets[0]["directory"], "Compare_transmissibilityW_vs_realtime_throat_46.png")
plt.savefig(plot_filename, bbox_inches='tight')
plt.show()
