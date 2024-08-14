import os
import numpy as np

in_folder_path = "data/K562"
for file_name in os.listdir(in_folder_path):
    file_path = os.path.join(in_folder_path, file_name)
    print(file_path)

    # read file
    data = np.loadtxt(file_path, delimiter=None)
    interaction_counts = np.sum(data > 1, axis=1)
    valid_interactions = data * (data > 1)
    # calculate HiC values
    interaction_sums = interaction_counts * np.sum(valid_interactions, axis=1)

    # save file
    output_file_path = f"prods/K562/{file_name}.txt"
    np.savetxt(output_file_path, interaction_sums, fmt="%.3f")