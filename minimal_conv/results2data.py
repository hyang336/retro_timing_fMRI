import os
import shutil
import csv

# Define subject list
sub_list='/home/users/hyang336/retro_timing_fMRI/sub_list_yaam'

# Define run list
run_list='/home/users/hyang336/retro_timing_fMRI/run_list'

# Define the source directory where the nifti files are located
source_dir = '/scratch/users/hyang336/retro_timing_fMRI_data/yaam_backup20240124'

# Define the destination directory where the files will be moved and renamed
destination_dir = '/scratch/users/hyang336/retro_timing_fMRI_data/yaam_torch_data'

# Create the destination directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# Initialize a list to store the annotations
annotations = []

# Read the subject and run lists
with open(sub_list, 'r') as file:
    sub_list = file.readlines()

with open(run_list, 'r') as file:
    run_list = file.readlines()

# Iterate over subjects, read lines from the subject list    
for participant in sub_list:
    participant = participant.strip()
    print(f'Processing participant {participant}')
    
    # Iterate over runs, read lines from the run list    
    for run in run_list:
        run = run.strip()
        print(f'Processing run {run}')

        # Check if the directory is not empty
        if os.listdir(f'{source_dir}/{participant}_Run_{run}_ResMS_fold-1'):
            # Iterate over the files in the directory
            for filename in os.listdir(f'{source_dir}/{participant}_Run_{run}_ResMS_fold-1'):
                # Check if the file is a nifti file with a given pattern
                if filename.startswith('ResMS') and filename.endswith('.nii'):
                    # Extract the time-shift value between ResMS and .nii from the filename
                    time_shift = filename.split('ResMS')[1].split('.nii')[0]
                    # convert the time-shift value to a float
                    time_shift = float(time_shift)

                    # Generate the new filename by affixating the participant and run information onto the original filename
                    new_filename = f'{participant}_run{run}_time{time_shift}.nii'

                    # Copy and rename the file to the destination directory
                    shutil.copyfile(f'{source_dir}/{participant}_Run_{run}_ResMS_fold-1/{filename}', os.path.join(destination_dir, new_filename))

                    # Append the annotation to the list
                    annotations.append([new_filename, time_shift])

# Write the annotations to a csv file
csv_file = f'{destination_dir}/annotations.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Filename', 'time_shift'])
    writer.writerows(annotations)

print('Data preparation completed.')
