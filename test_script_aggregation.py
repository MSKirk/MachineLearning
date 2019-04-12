""" Aggregation script.

Maintained at https://github.com/MSKirk/MachineLearning/blob/master/script_aggregation.py

Merges the content of the csv files of all YEAR_MONTH subdirectories into a single "global" csv file (global_csv_file)
and move all files to a common jp2 and labels subdirectory without duplicates.
User must update parent_dir, parent_dir2 and global_csv_file to personal case.
Note that depending on how you got the data, you may already have the global csv file under parent_dir.
It will be overwritten upon execution of this script, unless its path (global_csv_file) is changed.

Given a parent_directory hosting the original data (parent_dir):

parent_dir
    2010_12
        jp2
        labels_masks
    2011_01
        jp2
        labels_masks
    ...
    2019_01
        jp2
        labels_masks
    empty_sets (please ignore this)
    label_jp2_map_global.csv

Upon completion of this script we end up for a parent directory hosting the moved duplicate-free data (parent_dir2):

parent_dir2
    jp2
    labels_masks
    global_csv_file

Where global_csv_file is the aggregated csv file without duplicated file paths.

Note that this script is non-destructive and will only move data without deleting the duplicated files left in the
original directories.

In case no backup is made, 2 conveniency csv files are created to move the data back to the original directory tree,
if ever needed.

This script is of course optional as long as you can map properly map the data in all subdirectories in a way that accounts
for the duplicates across consecutive months.

Authors: Dr. Raphael Attie (raphael.attie@nasa.gov) & Dr. Michael Kirk (michael.s.kirk@nasa.gov)
"""

import os
import glob
import pandas as pd
import shutil
import csv


#############   Set some data directories - update to your personal case  #############

# Parent directory of all YEAR_MONTH subdirectories that will also contain the global csv file
parent_dir = '/Users/rattie/temp/'#'/Volumes/SolarData/LabeledImages/'
# Common directory where all files will be moved, without duplicates.
parent_dir2 = parent_dir
# Filename of csv file that will be the aggregation all csv files of all YEAR_MONTH subdirectories without duplicates
global_csv_file = os.path.join(parent_dir, 'label_jp2_map_global.csv')

######### (1) Creating the aggregated map of jp2 and label masks ###########

# Fetch the csv file paths recursively
#csv_files = sorted(glob.glob(os.path.join(parent_dir, '20*/label_jp2_map.csv')))
csv_files = sorted(glob.glob(os.path.join(parent_dir, 'temp*/label_jp2_map.csv')))
# Read their content and concatenate in a unique dataframe
dfs = []
for csvf in csv_files:
    print(csvf)
    dfs.append(pd.read_csv(csvf, header=None))

# Concatenate the dataframes into a single one while dropping all duplicates
label_jp2_map_global = pd.concat(dfs).drop_duplicates().reset_index(drop=True)
headers = ['npz file', 'HEK time', 'jp2 AIA 1600', 'jp2 AIA 1700', 'jp2 AIA 94', 'jp2 AIA 131', 'jp2 AIA 171',
           'jp2 AIA 193', 'jp2 AIA 211', 'jp2 AIA 304', 'jp2 AIA 335', 'jp2 HMI continuum', 'jp2 HMI magnetogram']

# Concatenate file basename with their parent directory relative path


# Write the new dataframe to a csv file
label_jp2_map_global.to_csv(global_csv_file, index=False, header=headers)



######### (2) Map the file paths for moving to a common directory #########
# Create csv of jp2 files to map the content of YEAR_MONTH-based paths to a common directory and vice versa,
# if moving them there is ever needed. Paths are written relative to parent directories.

# Setup common destination directory
jp2_dir = os.path.join(parent_dir2, 'jp2')
labels_dir = os.path.join(parent_dir2, 'label_masks')
png_dir = os.path.join(labels_dir, 'png')

# if not os.path.exists(jp2_dir):
#     os.makedirs(jp2_dir)
#
# if not os.path.exists(labels_dir):
#     os.makedirs(labels_dir)
#
# if not os.path.exists(png_dir):
#     os.makedirs(png_dir)

# Paths to original files
# jp2f = sorted(glob.glob(os.path.join(parent_dir2, '20*/jp2/*.jp2')))
# labels = sorted(glob.glob(os.path.join(parent_dir2, '20*/label_masks/*.*')))
jp2f = sorted(glob.glob(os.path.join(parent_dir, 'temp*/jp2/*.jp2')))
labels = sorted(glob.glob(os.path.join(parent_dir, 'temp*/label_masks/*.*')))

# List of original file paths free of duplicates to register in csv files
jp2f_list = []
labels_list = []
jp2f_csv = os.path.join(parent_dir2, 'map_non_duplicated_jp2_paths.csv')
labels_csv = os.path.join(parent_dir2, 'map_non_duplicated_labels_paths.csv')

# Map the jp2 files
new_files = []
for file in jp2f:
    new_file = os.path.join(jp2_dir, os.path.basename(file))
    print(new_file)
    if not new_file in new_files:
        # Map relative paths with respect to parent_directories
        original_file_relative = os.path.relpath(file, parent_dir)
        new_file_relative = os.path.relpath(new_file, parent_dir2)
        jp2f_list.append([original_file_relative, new_file_relative])
        #shutil.move(file, jp2_dir)
        new_files.append(new_file)


with open(jp2f_csv, 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(jp2f_list)
csvFile.close()

# Map the label masks (.npz files)
new_files = []
for file in labels:
    original_file_relative = os.path.relpath(file, parent_dir)
    _, ext = os.path.splitext(file)
    if ext == '.npz':
        new_file = os.path.join(labels_dir, os.path.basename(file))
    else:
        new_file = os.path.join(png_dir, os.path.basename(file))
    print(new_file)
    new_file_relative = os.path.relpath(new_file, parent_dir)
    if not new_file in new_files:
        labels_list.append([new_file, file])
        #shutil.move(file, new_file)
        new_files.append(new_file)

# Create the restore csv of .npz files (including png files) to move back to original paths if needed
with open(labels_csv, 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(labels_list)
csvFile.close()


