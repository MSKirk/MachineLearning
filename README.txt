To handle the duplicated/overlapping events and associated files accross consecutive months,
an aggregation script is available to map and move the data in a way that will get rid of these overlaps.
You can get get this script at https://github.com/MSKirk/MachineLearning/blob/master/script_aggregation.py
It can create (1) an aggregated map in a csv file and (2) move the data into a single parent directory while getting rid of the overlaps.

Here we only have performed (1): such aggregated map is available in the parent directory of all the YEAR_MONTH subidrectories.
However due to the great size of this data series and some instability of our usb drive connectivity that is not (yet) backed up,
we prefer to not perform (2) and leave the original tree as is instead of moving the nearly half TB of jp2, .npz, and png into a single "global" directory,
and rather let you decide how/where you want to move the files.