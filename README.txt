To handle the overlapping events and associated duplicated files accross consecutive months,
you will find:

(1) an aggregated map of all duplicate-free events in a csv file. This only shows file basenames:
label_jp2_map_global.csv

(2) csv files for mappin, copying or moving duplicate-free data relative paths into a single parent directory:
- map_non_duplicated_jp2_paths.csv
- map_non_duplicated_labels_paths.csv


The aggregation script to produce these csv files using the data on the original disk is available
at https://github.com/MSKirk/MachineLearning/blob/master/script_aggregation.py


