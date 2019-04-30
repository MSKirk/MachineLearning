To handle the overlapping events and associated duplicated files accross consecutive months,
you will find:
(1) an aggregated map of all overlap-free events in a csv file. This only shows file basenames
(2) csv files for mapping the data relative paths into a single parent directory.
All these maps get rid of the overlapping events and duplicated files by keeping only one
wherever they occur.

The aggregation script to produce these csv files is available
at https://github.com/MSKirk/MachineLearning/blob/master/script_aggregation.py


