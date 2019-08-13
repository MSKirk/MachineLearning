import os
from read_jp2 import read_solar_jp2
import matplotlib
matplotlib.use('Tkagg')
import matplotlib.pyplot as plt
import numpy as np
import calibration

# Get a jp2 sample. here using the jp2 included in this github repo
filepath = '../images/2011_06_25__00_59_43_71__SDO_AIA_AIA_1700.jp2'

pdata, pheader = read_solar_jp2(filepath)
pdata[0:1000, :] = 0

rdata = calibration.scale_rotate(pdata, 45)

# Display the image and make sure it's in the correct orientation with respect to the png sample
vmax = np.percentile(pdata, 99.5)

fs = 20
plt.figure(0, figsize=(18, 18))
plt.subplot(2, 2, 1)
plt.imshow(pdata, vmin=pdata.min(), vmax=vmax, origin='lower', cmap='gray')
plt.title('origin lower, no rotation', fontsize=fs)
plt.subplot(2, 2, 2)
plt.imshow(rdata, vmin=pdata.min(), vmax=vmax, origin='lower', cmap='gray')
plt.title('origin lower, rotation argument +45 deg', fontsize=fs)
plt.subplot(2, 2, 3)
plt.imshow(pdata, vmin=pdata.min(), vmax=vmax, cmap='gray')
plt.title('origin top, no rotation', fontsize=fs)
plt.subplot(2, 2, 4)
plt.imshow(rdata, vmin=pdata.min(), vmax=vmax, cmap='gray')
plt.title('origin top, rotation argument +45 deg', fontsize=fs)
plt.tight_layout()
plt.show()