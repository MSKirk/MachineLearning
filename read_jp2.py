import numpy as np
from astropy import time
from astropy import units as u
import sunpy.io
import calibration
import warnings

# initialize effective area class to avoid rereading the calibration table

aia_effective_area = calibration.AIAEffectiveArea()


def read_solar_jp2(filepath, verbose=False):

    img = sunpy.io.read_file(filepath, filetype='jp2')[0]

    if img.header['CROTA2'] != 0:
        if verbose:
            print('Rotating image to solar north')
        prepped_data = calibration.scale_rotate(img.data, img.header['CROTA2'])
    else:
        prepped_data = img.data

    if verbose:
        print('Correcting for CCD degradation')

    prepped_data = prepped_data * (1./aia_effective_area.effective_area_ratio(img.header['WAVELNTH']*u.AA, time.Time(img.header['DATE-OBS']).to_datetime()))

    if img.header['WAVELNTH'] in [94, 131, 171, 193, 211, 304, 335]:
        if np.sum(img.data == 0) / np.sum(img.data != 0) > 0.5:
            warnings.warn('Significant amount of data missing in the image', UserWarning)
    else:
        if np.sum(img.data == 0) / np.sum(img.data != 0) > 1.2:
            warnings.warn('Significant amount of data missing in the image', UserWarning)

    return prepped_data




