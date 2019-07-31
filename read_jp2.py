import numpy as np
from astropy import time
from astropy import units as u
import sunpy.io
import calibration
import warnings

# initialize effective area class to avoid rereading the calibration table
aia_effective_area = calibration.AIAEffectiveArea()

# The aia image size is fixed by the size of the detector. For AIA raw data, this has no reason to change.
aia_image_size = 4096


def read_solar_jp2(filepath, verbose=False):
    '''
    :param filepath: The full file path of the SDO .jp2 image to be read in
    :param verbose: Boolean, if True will print status statements
    :return: numpy array of prepped image
    '''

    # Read the image and header
    img = sunpy.io.read_file(filepath, filetype='jp2')[0]
    prepped_header = img.header

    # Rotation of image to Solar North
    if img.header['CROTA2'] != 0:
        if verbose:
            print('Rotating image to solar north')
        prepped_data = calibration.scale_rotate(img.data, img.header['CROTA2'])

        center = ((np.array(prepped_data.shape) - 1) / 2.0).astype(int)
        half_size = int(aia_image_size / 2)
        prepped_data = prepped_data[center[1] - half_size:center[1] + half_size, center[0] - half_size:center[0] + half_size]
        prepped_header['CROTA2'] = 0
    else:
        prepped_data = img.data

    # Normalizing the image intensity to levels at the start of the mission for AIA
    if 'AIA' in img.header['INSTRUME']:
        if verbose:
            print('Correcting for CCD degradation')

        prepped_data *= (1./aia_effective_area.effective_area_ratio(img.header['WAVELNTH']*u.AA, time.Time(img.header['DATE-OBS']).to_datetime()))
        prepped_data[prepped_data < 0] = 0
        prepped_header['DATAMIN'] = 0

    # User Warning if there is significant amount of missing data
    non_zero = np.count_nonzero(prepped_data != 0)

    if img.header['WAVELNTH'] in [94, 131, 171, 193, 211, 304, 335]:
        if (prepped_data.size - non_zero) / non_zero > 0.6:
            warnings.warn('Significant amount of data missing in the image', UserWarning)
    else:
        if (prepped_data.size - non_zero) / non_zero > 1.2:
            warnings.warn('Significant amount of data missing in the image', UserWarning)

    return prepped_data, prepped_header

def write_solar_jp2(fname, image_array, image_header):

    sunpy.io.write_file(fname, image_array, image_header, filetype='auto')


