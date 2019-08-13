import os
from read_jp2 import read_solar_jp2
import matplotlib
matplotlib.use('Tkagg')

# Get a jp2 sample. here using the jp2 included in this github repo
filepath = '../images/2011_06_25__00_59_43_71__SDO_AIA_AIA_1700.jp2'

pdata = read_solar_jp2(filepath)
