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

#Extract p_gas data from a VTP file
def extract_p_gas_from_vtp(file_path):
    try:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        polydata = reader.GetOutput()
        p_gas_array = np.array(polydata.GetPointData().GetArray("p_gas"))
        return p_gas_array
    except Exception as e:
        print(f"Failed to extract p_gas data from {file_path}: {e}")
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
    
    p_gas_data_arrays = [extract_p_gas_from_vtp(file_path) for file_path in file_paths]
    p_gas_data_arrays = [data for data in p_gas_data_arrays if data is not None]
    
    return timesteps, p_gas_data_arrays

# Directories and file patterns for each dataset
datasets = [
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_sink/time_1sec/", "pvd": "2p_im_magic_time_1sec.pvd", "base_filename": "2p_im_magic_time_1sec"},
    {"directory": "/home/hp/dumux/harsha_pnm/build-cmake/appl/2p_magic/time_1sec/", "pvd": "2p_magic_reg_time_1sec.pvd", "base_filename": "2p_magic_reg_time_1sec"}
]

# Initialize lists to hold the merged data
merged_timesteps = []
merged_p_gas_data = []
labels = []

# Loop through each dataset and extract data
for dataset in datasets:
    timesteps, p_gas_data_arrays = extract_data_from_dataset(
        dataset["directory"], dataset["pvd"], dataset["base_filename"])
    
    if timesteps and p_gas_data_arrays:
        merged_timesteps.append(timesteps)
        merged_p_gas_data.append(p_gas_data_arrays)
        labels.append(dataset["base_filename"])
    else:
        print(f"Failed to process dataset: {dataset['pvd']}")

# Plotting the merged data for pore 
pore_index = 22
plt.figure(figsize=(10, 6))
plt.xlabel('Time [s]', fontsize=12)
plt.ylabel('Local $p_{n}$ [Pa]', fontsize=12)
plt.title(f'Pore {pore_index} Comparison', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=10)
plt.tick_params(axis='both', which='minor', direction='out', length=3, width=1)
plt.gca().xaxis.set_tick_params(width=1.5)
plt.gca().yaxis.set_tick_params(width=1.5)


colors = ['blue', 'red']
markers = ['o','*']

for i, (timesteps, p_gas_data_arrays) in enumerate(zip(merged_timesteps, merged_p_gas_data)):
    for j, p_gas_array in enumerate(p_gas_data_arrays):
        time_step = timesteps[j]
        plt.plot(time_step, p_gas_array[pore_index], color=colors[i], marker=markers[i], markersize=6, label=f'{labels[i]}' if j == 0 else "")

plt.legend()
plt.tight_layout()

# Save the plot as a .png file
plot_filename = os.path.join(datasets[0]["directory"], "Compare_p_gas_vs_realtime_pore_22.png")
plt.savefig(plot_filename, bbox_inches='tight')
plt.show()
