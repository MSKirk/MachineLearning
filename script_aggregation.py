import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY
import datetime
from sunpy.time import parse_time
import glob
import pandas as pd
import numpy as np


parent_dir = '/Volumes/SolarData/LabledImages/'
# Fetch the csv file paths recursively
csv_files = sorted(glob.glob(os.path.join(parent_dir, '20*/label_jp2_map.csv')))
# Read their content and concatenate in a unique dataframe
#df = pd.read_csv(csv_files[0])
dfs = []
for csvf in csv_files:
    print(csvf)
    dfs.append(pd.read_csv(csvf, header=None))

label_jp2_map_global = pd.concat(dfs).drop_duplicates().reset_index(drop=True)
headers = ['npz file', 'HEK time', 'jp2 AIA 1600', 'jp2 AIA 1700', 'jp2 AIA 94', 'jp2 AIA 131', 'jp2 AIA 171',
           'jp2 AIA 193', 'jp2 AIA 211', 'jp2 AIA 304', 'jp2 AIA 335', 'jp2 HMI continuum', 'jp2 HMI magnetogram']

label_jp2_map_global.to_csv(os.path.join(parent_dir, 'label_jp2_map_global.csv'), index=False, header=headers)
