import numpy as np
from astropy import time
from astropy import units as u
import sunpy.io
import calibration

# pull out effective area class to avoid rereading the table

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

    prepped_data /= aia_effective_area.effective_area_ratio(img.header['WAVELNTH']*u.AA, time.Time(img.header['DATE-OBS']).to_datetime())



    return prepped_data




